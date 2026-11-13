/**
 * Draws the chamber with WebGL2. Every lit pixel is a droplet: each is a point whose brightness
 * is its weight times its cross-section (r^2) times how strongly a droplet scatters light from
 * the lamp toward the camera.
 *
 * Droplets are tens of wavelengths across, so they scatter mostly forward. That is modelled with
 * a Henyey-Greenstein phase function with asymmetry g = 0.85, the usual value for cloud droplets
 * in visible light; it is why tracks look brightest when you look toward the lamp.
 */
import type { Engine } from './engine'

export interface Camera {
  /** 'orbit' is a perspective view; 'top' looks straight down with no perspective, for measuring. */
  mode: 'orbit' | 'top'
  /** Orbit: angle around the vertical axis and elevation above the plate, radians. */
  yaw: number
  pitch: number
  /** Orbit: distance from the target, m. */
  distance: number
  target: [number, number, number]
}

export interface DrawOptions {
  /** Brightness multiplier. */
  exposure: number
  /** Sprite diameter in device pixels. */
  sprite: number
  /** Tint droplets by the particle that made them (not physical; droplets are all white). */
  colourByParticle: boolean
  /** Draw the plate and the chamber edges. */
  scenery: boolean
}

export const DEFAULT_DRAW: DrawOptions = {
  exposure: 1,
  sprite: 3,
  colourByParticle: false,
  scenery: true,
}

export const ASYMMETRY = 0.85

const DROPLET_VS = `#version 300 es
precision highp float;
layout(location = 0) in vec3 a_pos;
layout(location = 1) in float a_radius;
layout(location = 2) in float a_weight;
layout(location = 3) in float a_kind;
uniform mat4 u_viewProj;
uniform vec3 u_eye;
uniform vec3 u_eyeDir;   // used instead of u_eye when u_ortho = 1
uniform float u_ortho;
uniform vec3 u_light;    // direction the light travels
uniform float u_exposure;
uniform float u_sprite;
uniform float u_tint;
out vec3 v_colour;
const float G = ${ASYMMETRY.toFixed(2)};
float phase(float c) { return (1.0 - G * G) / pow(1.0 + G * G - 2.0 * G * c, 1.5); }
void main() {
  if (a_radius <= 0.0) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); gl_PointSize = 1.0; v_colour = vec3(0.0); return; }
  gl_Position = u_viewProj * vec4(a_pos, 1.0);
  gl_PointSize = u_sprite;
  vec3 toEye = mix(normalize(u_eye - a_pos), -u_eyeDir, u_ortho);
  // Brightness relative to one 10-micron droplet seen at right angles to the beam.
  float r = a_radius / 1.0e-5;
  float b = u_exposure * a_weight * r * r * phase(dot(u_light, toEye)) / phase(0.0);
  b *= 9.0 / (u_sprite * u_sprite); // same total light whatever the sprite size
  vec3 white = vec3(0.86, 0.93, 1.0);
  vec3 tint = a_kind < 0.5 ? vec3(0.55, 0.65, 0.8) : a_kind < 1.5 ? vec3(1.0, 0.62, 0.3) : a_kind < 2.5 ? vec3(0.45, 1.0, 0.6) : vec3(0.6, 0.75, 1.0);
  v_colour = b * mix(white, tint, u_tint);
}`

const DROPLET_FS = `#version 300 es
precision highp float;
in vec3 v_colour;
out vec4 o;
void main() {
  vec2 d = gl_PointCoord * 2.0 - 1.0;
  float fall = exp(-3.0 * dot(d, d));
  o = vec4(v_colour * fall, 1.0);
}`

const FLAT_VS = `#version 300 es
precision highp float;
layout(location = 0) in vec3 a_pos;
uniform mat4 u_viewProj;
out vec3 v_world;
void main() { v_world = a_pos; gl_Position = u_viewProj * vec4(a_pos, 1.0); }`

const FLAT_FS = `#version 300 es
precision highp float;
in vec3 v_world;
uniform vec4 u_colour;
uniform float u_plate;
out vec4 o;
void main() {
  // The plate is black anodized metal: almost nothing, a little brighter toward the lamp.
  float sheen = u_plate * 0.02 * smoothstep(-0.12, 0.12, v_world.x);
  o = vec4(u_colour.rgb + sheen, u_colour.a);
}`

function compile(gl: WebGL2RenderingContext, type: number, source: string): WebGLShader {
  const shader = gl.createShader(type)
  if (!shader) throw new Error('cannot create shader')
  gl.shaderSource(shader, source)
  gl.compileShader(shader)
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS))
    throw new Error(`shader compile: ${gl.getShaderInfoLog(shader)}`)
  return shader
}

function link(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
  const program = gl.createProgram()
  if (!program) throw new Error('cannot create program')
  gl.attachShader(program, compile(gl, gl.VERTEX_SHADER, vs))
  gl.attachShader(program, compile(gl, gl.FRAGMENT_SHADER, fs))
  gl.linkProgram(program)
  if (!gl.getProgramParameter(program, gl.LINK_STATUS))
    throw new Error(`program link: ${gl.getProgramInfoLog(program)}`)
  return program
}

type Mat4 = Float32Array
type Vec3 = [number, number, number]

function multiply(a: Mat4, b: Mat4): Mat4 {
  const out = new Float32Array(16)
  for (let c = 0; c < 4; c++)
    for (let r = 0; r < 4; r++) {
      let sum = 0
      for (let k = 0; k < 4; k++) sum += a[k * 4 + r] * b[c * 4 + k]
      out[c * 4 + r] = sum
    }
  return out
}

function lookAt(eye: Vec3, target: Vec3, up: Vec3): Mat4 {
  const f = normalize([target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]])
  const s = normalize(cross(f, up))
  const u = cross(s, f)
  // biome-ignore format: matrix layout
  return new Float32Array([
    s[0], u[0], -f[0], 0,
    s[1], u[1], -f[1], 0,
    s[2], u[2], -f[2], 0,
    -dot(s, eye), -dot(u, eye), dot(f, eye), 1,
  ])
}

function perspective(fovY: number, aspect: number, near: number, far: number): Mat4 {
  const f = 1 / Math.tan(fovY / 2)
  // biome-ignore format: matrix layout
  return new Float32Array([
    f / aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, (far + near) / (near - far), -1,
    0, 0, (2 * far * near) / (near - far), 0,
  ])
}

function orthographic(halfW: number, halfH: number, near: number, far: number): Mat4 {
  // biome-ignore format: matrix layout
  return new Float32Array([
    1 / halfW, 0, 0, 0,
    0, 1 / halfH, 0, 0,
    0, 0, -2 / (far - near), 0,
    0, 0, -(far + near) / (far - near), 1,
  ])
}

const dot = (a: Vec3, b: Vec3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
const cross = (a: Vec3, b: Vec3): Vec3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
]
function normalize(v: Vec3): Vec3 {
  const n = Math.hypot(v[0], v[1], v[2])
  return [v[0] / n, v[1] / n, v[2] / n]
}

/** Margin around the box in the top view, as a fraction of its half-width. */
export const TOP_VIEW_MARGIN = 1.05

/** Pixels per metre in the top view, for a canvas `width` pixels wide over a box of half-width `halfX`. */
export function topViewScale(width: number, halfX: number): number {
  return width / (2 * halfX * TOP_VIEW_MARGIN)
}

export class Renderer {
  readonly gl: WebGL2RenderingContext
  private readonly droplets: WebGLProgram
  private readonly flat: WebGLProgram
  private readonly vao: WebGLVertexArrayObject
  private readonly buffers: {
    pos: WebGLBuffer
    radius: WebGLBuffer
    weight: WebGLBuffer
    kind: WebGLBuffer
  }
  private readonly sceneryVao: WebGLVertexArrayObject
  private readonly sceneryBuffer: WebGLBuffer
  private capacity = 0
  /** Direction the lamp's light travels: from the +x wall, level with the plate. */
  readonly light: Vec3 = normalize([-1, 0.15, -0.05])

  constructor(canvas: HTMLCanvasElement | OffscreenCanvas) {
    const gl = canvas.getContext('webgl2', {
      alpha: false,
      antialias: false,
      preserveDrawingBuffer: true,
    }) as WebGL2RenderingContext | null
    if (!gl) throw new Error('WebGL2 is not available')
    this.gl = gl
    this.droplets = link(gl, DROPLET_VS, DROPLET_FS)
    this.flat = link(gl, FLAT_VS, FLAT_FS)
    const make = () => {
      const b = gl.createBuffer()
      if (!b) throw new Error('cannot create buffer')
      return b
    }
    this.buffers = { pos: make(), radius: make(), weight: make(), kind: make() }
    const vao = gl.createVertexArray()
    const sceneryVao = gl.createVertexArray()
    if (!vao || !sceneryVao) throw new Error('cannot create vertex array')
    this.vao = vao
    this.sceneryVao = sceneryVao
    gl.bindVertexArray(vao)
    const attrib = (buffer: WebGLBuffer, location: number, size: number, type: number) => {
      gl.bindBuffer(gl.ARRAY_BUFFER, buffer)
      gl.enableVertexAttribArray(location)
      gl.vertexAttribPointer(location, size, type, false, 0, 0)
    }
    attrib(this.buffers.pos, 0, 3, gl.FLOAT)
    attrib(this.buffers.radius, 1, 1, gl.FLOAT)
    attrib(this.buffers.weight, 2, 1, gl.FLOAT)
    attrib(this.buffers.kind, 3, 1, gl.UNSIGNED_BYTE)
    this.sceneryBuffer = make()
    gl.bindVertexArray(sceneryVao)
    attrib(this.sceneryBuffer, 0, 3, gl.FLOAT)
    gl.bindVertexArray(null)
  }

  /** View-projection matrix and eye for a camera and viewport aspect. */
  matrices(camera: Camera, aspect: number, halfX: number, halfY: number, height: number) {
    if (camera.mode === 'top') {
      const halfW = halfX * TOP_VIEW_MARGIN
      const eye: Vec3 = [0, 0, height + 1]
      const view = lookAt(eye, [0, 0, 0], [0, 1, 0])
      return {
        viewProj: multiply(orthographic(halfW, halfW / aspect, 0.1, height + 2), view),
        eye,
        eyeDir: [0, 0, -1] as Vec3,
        ortho: 1,
      }
    }
    const cp = Math.cos(camera.pitch)
    const eye: Vec3 = [
      camera.target[0] + camera.distance * cp * Math.cos(camera.yaw),
      camera.target[1] + camera.distance * cp * Math.sin(camera.yaw),
      camera.target[2] + camera.distance * Math.sin(camera.pitch),
    ]
    const view = lookAt(eye, camera.target, [0, 0, 1])
    const reach = camera.distance + Math.hypot(halfX, halfY, height)
    return {
      viewProj: multiply(perspective(0.6, aspect, 0.02, reach * 2), view),
      eye,
      eyeDir: [0, 0, -1] as Vec3,
      ortho: 0,
    }
  }

  draw(engine: Engine, camera: Camera, options: DrawOptions = DEFAULT_DRAW): void {
    const { gl } = this
    const width = gl.drawingBufferWidth
    const height = gl.drawingBufferHeight
    const { halfX, halfY, height: boxHeight } = engine.box
    const m = this.matrices(camera, width / height, halfX, halfY, boxHeight)
    gl.viewport(0, 0, width, height)
    gl.clearColor(0, 0, 0, 1)
    gl.clear(gl.COLOR_BUFFER_BIT)
    gl.disable(gl.DEPTH_TEST)

    if (options.scenery) this.drawScenery(m.viewProj, halfX, halfY, boxHeight)

    const n = engine.count
    gl.bindVertexArray(this.vao)
    const upload = (buffer: WebGLBuffer, data: ArrayBufferView, perItem: number, bytes: number) => {
      gl.bindBuffer(gl.ARRAY_BUFFER, buffer)
      if (this.capacity < engine.options.capacity)
        gl.bufferData(gl.ARRAY_BUFFER, engine.options.capacity * perItem * bytes, gl.DYNAMIC_DRAW)
      gl.bufferSubData(gl.ARRAY_BUFFER, 0, data, 0, n * perItem)
    }
    upload(this.buffers.pos, engine.pos, 3, 4)
    upload(this.buffers.radius, engine.radius, 1, 4)
    upload(this.buffers.weight, engine.weight, 1, 4)
    upload(this.buffers.kind, engine.kind, 1, 1)
    this.capacity = engine.options.capacity

    gl.useProgram(this.droplets)
    const u = (name: string) => gl.getUniformLocation(this.droplets, name)
    gl.uniformMatrix4fv(u('u_viewProj'), false, m.viewProj)
    gl.uniform3fv(u('u_eye'), m.eye)
    gl.uniform3fv(u('u_eyeDir'), m.eyeDir)
    gl.uniform1f(u('u_ortho'), m.ortho)
    gl.uniform3fv(u('u_light'), this.light)
    gl.uniform1f(u('u_exposure'), options.exposure)
    gl.uniform1f(u('u_sprite'), options.sprite)
    gl.uniform1f(u('u_tint'), options.colourByParticle ? 1 : 0)
    gl.enable(gl.BLEND)
    gl.blendFunc(gl.ONE, gl.ONE)
    if (n > 0) gl.drawArrays(gl.POINTS, 0, n)
    gl.disable(gl.BLEND)
    gl.bindVertexArray(null)
  }

  private drawScenery(viewProj: Mat4, hx: number, hy: number, h: number): void {
    const { gl } = this
    // biome-ignore format: vertex lists
    const plate = [-hx, -hy, 0, hx, -hy, 0, hx, hy, 0, -hx, -hy, 0, hx, hy, 0, -hx, hy, 0]
    const corners: Vec3[] = [
      [-hx, -hy, 0],
      [hx, -hy, 0],
      [hx, hy, 0],
      [-hx, hy, 0],
    ]
    const edges: number[] = []
    for (let i = 0; i < 4; i++) {
      const a = corners[i]
      const b = corners[(i + 1) % 4]
      edges.push(...a, ...b, a[0], a[1], h, b[0], b[1], h, ...a, a[0], a[1], h)
    }
    gl.bindVertexArray(this.sceneryVao)
    gl.bindBuffer(gl.ARRAY_BUFFER, this.sceneryBuffer)
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([...plate, ...edges]), gl.DYNAMIC_DRAW)
    gl.useProgram(this.flat)
    gl.uniformMatrix4fv(gl.getUniformLocation(this.flat, 'u_viewProj'), false, viewProj)
    gl.uniform4f(gl.getUniformLocation(this.flat, 'u_colour'), 0.012, 0.014, 0.018, 1)
    gl.uniform1f(gl.getUniformLocation(this.flat, 'u_plate'), 1)
    gl.drawArrays(gl.TRIANGLES, 0, 6)
    gl.uniform4f(gl.getUniformLocation(this.flat, 'u_colour'), 0.09, 0.1, 0.12, 1)
    gl.uniform1f(gl.getUniformLocation(this.flat, 'u_plate'), 0)
    gl.drawArrays(gl.LINES, 6, edges.length / 3)
    gl.bindVertexArray(null)
  }

  /** Read the frame back as RGBA bytes, top row first. */
  readPixels(): { width: number; height: number; data: Uint8Array } {
    const { gl } = this
    const width = gl.drawingBufferWidth
    const height = gl.drawingBufferHeight
    const raw = new Uint8Array(width * height * 4)
    gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, raw)
    const data = new Uint8Array(raw.length)
    for (let y = 0; y < height; y++)
      data.set(raw.subarray((height - 1 - y) * width * 4, (height - y) * width * 4), y * width * 4)
    return { width, height, data }
  }
}

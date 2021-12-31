const $ = require('jquery');
const dat = require('dat.gui');
import GPU from 'gpu.js';
import SimplexNoise from 'simplex-noise';
const seedrandom = require('seedrandom');


const addNoiseFunctions = function(GPU) {
	//
	// Description : Array and textureless GLSL 2D/3D/4D simplex
	//               noise functions.
	//      Author : Ian McEwan, Ashima Arts.
	//  Maintainer : ijm
	//     Lastmod : 20110822 (ijm)
	//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
	//               Distributed under the MIT License. See LICENSE file.
	//               https://github.com/ashima/webgl-noise
	//

//  GPU.addNativeFunction('mod2893',
//  `
//    vec3 mod2893(vec3 x) {
//      return x - floor(x * (1.0 / 289.0)) * 289.0;
//    }
//  `
//  );

//  GPU.addNativeFunction('mod2894',
//  `
//  vec4 mod2894(vec4 x) {
//    return x - floor(x * (1.0 / 289.0)) * 289.0;
//  }
//  `);

//  GPU.addNativeFunction('permute',
//  `
//  vec4 permute(vec4 x) {
//       return mod2894(((x*34.0)+1.0)*x);
//  }
//  `);

//  GPU.addNativeFunction('taylorInvSqrt',
//  `
//  vec4 taylorInvSqrt(vec4 r)
//  {
//    return 1.79284291400159 - 0.85373472095314 * r;
//  }
//  `);

	GPU.addNativeFunction('snoise',
	`
	float snoise4(vec3 v) {
    return 0.0;
  }
	vec3 mod289(vec3 x) {
		return x - floor(x * (1.0 / 289.0)) * 289.0;
	}
	vec4 mod289(vec4 x) {
		return x - floor(x * (1.0 / 289.0)) * 289.0;
	}
	vec4 permute(vec4 x) {
			 return mod289(((x*34.0)+1.0)*x);
	}
	vec4 taylorInvSqrt(vec4 r)
	{
		return 1.79284291400159 - 0.85373472095314 * r;
	}
	float snoise(vec3 v)
		{
		const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
		const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

	// First corner
		vec3 i  = floor(v + dot(v, C.yyy) );
		vec3 x0 =   v - i + dot(i, C.xxx) ;

	// Other corners
		vec3 g = step(x0.yzx, x0.xyz);
		vec3 l = 1.0 - g;
		vec3 i1 = min( g.xyz, l.zxy );
		vec3 i2 = max( g.xyz, l.zxy );

		//   x0 = x0 - 0.0 + 0.0 * C.xxx;
		//   x1 = x0 - i1  + 1.0 * C.xxx;
		//   x2 = x0 - i2  + 2.0 * C.xxx;
		//   x3 = x0 - 1.0 + 3.0 * C.xxx;
		vec3 x1 = x0 - i1 + C.xxx;
		vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
		vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

	// Permutations
		i = mod289(i);
		vec4 p = permute( permute( permute(
							 i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
						 + i.y + vec4(0.0, i1.y, i2.y, 1.0 ))
						 + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

	// Gradients: 7x7 points over a square, mapped onto an octahedron.
	// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
		float n_ = 0.142857142857; // 1.0/7.0
		vec3  ns = n_ * D.wyz - D.xzx;

		vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

		vec4 x_ = floor(j * ns.z);
		vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

		vec4 x = x_ *ns.x + ns.yyyy;
		vec4 y = y_ *ns.x + ns.yyyy;
		vec4 h = 1.0 - abs(x) - abs(y);

		vec4 b0 = vec4( x.xy, y.xy );
		vec4 b1 = vec4( x.zw, y.zw );

		//vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
		//vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
		vec4 s0 = floor(b0)*2.0 + 1.0;
		vec4 s1 = floor(b1)*2.0 + 1.0;
		vec4 sh = -step(h, vec4(0.0));

		vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
		vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

		vec3 p0 = vec3(a0.xy,h.x);
		vec3 p1 = vec3(a0.zw,h.y);
		vec3 p2 = vec3(a1.xy,h.z);
		vec3 p3 = vec3(a1.zw,h.w);

	//Normalise gradients
		vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
		p0 *= norm.x;
		p1 *= norm.y;
		p2 *= norm.z;
		p3 *= norm.w;

	// Mix final noise value
		vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
		m = m * m;
		return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1),
																	dot(p2,x2), dot(p3,x3) ) );
	}
	`);
};

const createShape = function (simplex, ws, hs) {
  const colors = new Float32Array(4 * ws * hs);
  const seed = new Date().getTime();
  const initradius = Math.random() * 0.5;
  const spikiness = Math.random() * 1.0;
  for (let x = 0; x < ws; x++) {
    for (let y = 0; y < hs; y++) {
      const dx = x - 0.5 * ws;
      const dy = y - 0.5 * hs;
      const radius = 0.5 * Math.max(ws, hs);
      const vradius =
        (0.75 + 0.25 * Math.sin(10.0 * Math.atan2(dy, dx))) * radius;
      const lerp = (x, y, a) => x * (1 - a) + y * a;
      const dlerp = (x, a) => lerp(x, x, a);
      const angle = Math.atan2(dy, dx) + Math.PI;
      const nangle = angle / Math.PI / 2.0;
      const c = 1.0;
      const vradius2 = lerp(
        (initradius + spikiness * Math.abs(simplex.noise2D(c * angle, seed))) * radius,
        (initradius + spikiness * Math.abs(simplex.noise2D(c * -angle, seed))) * radius,
        nangle
      );
      const currentRadius = Math.sqrt(dx * dx + dy * dy);
      const innerShape = currentRadius <= vradius2;
      const border =
        currentRadius >= vradius2 && currentRadius <= vradius2 + 1.0;
      // const idx = 4 * (y * w + x);
      const bc = () =>
        0.5 * (simplex.noise3D((0.5 * x) / ws, (0.5 * y) / ws, seed) + 1.0);

      const r = innerShape ? 0.0 : 1.0;
      const g = border ? 1.0 : 0.0;
      const b = innerShape ? bc() : 0.0;
      const a = innerShape || border ? 1.0 : 0.0;
      const t = 4 * (y * ws + x);
      const shape = innerShape ? bc() : 0.0;
      colors[t + 0] = shape;
      colors[t + 1] = shape;
      colors[t + 2] = shape;
      colors[t + 3] = 1.0;
    }
  }
  return colors;
};

const main = function(options) {

  seedrandom(options.seed, { global: true });
  const simplex = new SimplexNoise('seed');

  console.log(options);
  const w = Math.floor(parseInt(options.width));
  const h = Math.floor(parseInt(options.height));

  const first = $("#mainImage").length == 0;

  const canvas = $("#mainImage")[0] || $('<canvas id="mainImage" />')[0];

  canvas.width = w;
  canvas.height = h;

  if (first) {
    const div = $("#main")[0] || $('<div id="main" />')[0];
    $("body")
      .css('margin', '0px')
      .css('padding', '0px')
      .append(div);
    $("#main")
      .css('width', '100%')
      .css('height', '100vh')
      .css('background', 'black')
      .css('display', 'flex')
      .css('justify-content', 'center')
      .css('align-items', 'center')
      .append(canvas);
  }

  const ctx = canvas.getContext('2d');
  $("#mainImage")
    .css('border', '1px solid white')

  const data = ctx.createImageData(w, h);

  const points = [];
  const pointCount = options.pointsEnabled ? options.pointCount : 1;
  const pointsX = new Float32Array(pointCount);
  const pointsY = new Float32Array(pointCount);
  const pointsColorR = new Float32Array(pointCount);
  const pointsColorG = new Float32Array(pointCount);
  const pointsColorB = new Float32Array(pointCount);

  for (let i = 0; i < pointCount; i++) {
    const c = {
      r: options.pointsRandomColor ? Math.random() : options.pointsColor[0],
      g: options.pointsRandomColor ? Math.random() : options.pointsColor[1],
      b: options.pointsRandomColor ? Math.random() : options.pointsColor[2],
    };
    const point = {
      x: Math.floor(Math.random() * w),
      y: Math.floor(Math.random() * h),
      color: c,
    };
    points.push(point);
    pointsX[i] = point.x;
    pointsY[i] = point.y;
    pointsColorR[i] = point.color.r;
    pointsColorG[i] = point.color.g;
    pointsColorB[i] = point.color.b;
  }

	const ws = options.shapeWidth;
	const hs = options.shapeHeight;
  const shapeCount = 100;
  const shapeLocations = new Float32Array(shapeCount * 2);
  const shapeColors = new Float32Array(shapeCount * 3);
	const shapes = new Float32Array(shapeCount * 4 * ws * hs);
  for (let i = 0; i < shapeCount; i++) {
    const shape = createShape(simplex, ws, hs);
    const color = {
      r: Math.random(),
      g: Math.random(),
      b: Math.random()
    };
    shapeLocations[i * 2] = Math.floor(Math.random() * w);
    shapeLocations[i * 2 + 1] = Math.floor(Math.random() * h);
    shapeColors[i * 3 + 0] = color.r;
    shapeColors[i * 3 + 1] = color.g;
    shapeColors[i * 3 + 2] = color.b;
    for (let x = 0; x < ws; x++) {
      for (let y = 0; y < hs; y++) {
        const idx = 4 * (y * ws + x);
        const gidx = 4 * i * ws * hs + idx;
        shapes[gidx + 0] = shape[idx + 0] * color.r;
        shapes[gidx + 1] = shape[idx + 1] * color.g;
        shapes[gidx + 2] = shape[idx + 2] * color.b;
        shapes[gidx + 3] = shape[idx + 3];
      }
    }
	}


  const gpu = new GPU.GPU();
	addNoiseFunctions(gpu);
  function lerp2(x, y, a) {
    return x * (1.0 - a) + y * a;
  }
  function snoise3(a, b, c) {
    const ret = snoise([a, b, c]);
    return ret;
  }
  function createShapeGPU(initRadius, spikiness, shapeSeed, ws, hs, x, y) {
    const dx = x - 0.5 * ws;
    const dy = y - 0.5 * hs;
    const radius = 0.5 * Math.max(ws, hs);
    const vradius =
      (0.75 + 0.25 * Math.sin(10.0 * Math.atan2(dy, dx))) * radius;
    const angle = Math.atan2(dy, dx) + Math.PI;
    const nangle = angle / Math.PI / 2.0;
    const c = 1.0;
    const from = (initradius + spikiness * Math.abs(snoise3(c * angle, shapeSeed, 0.0))) * radius;
    const to = (initradius + spikiness * Math.abs(snoise3(c * -angle, shapeSeed, 0.0))) * radius;
    const vradius2 = lerp2(from, to, nangle);
    const currentRadius = Math.sqrt(dx * dx + dy * dy);
    const border = currentRadius >= vradius2 && currentRadius <= vradius2 + 1.0;
    const bc1 = snoise3((0.5 * x) / ws, (0.5 * y) / ws, shapeSeed);
    const bc2 = 0.5 * bc1 + 1.0;
    const shape = currentRadius <= vradius2 ? bc2 : 0.0;
    return shape;
  };
	gpu.setFunctions([lerp2, snoise3, createShapeGPU]);

  const maxd = Math.sqrt(w * w + h * h);
  const maxe = 1.0 / pointCount;
  const process = gpu
    .createKernel(function (
      pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB,
      shapes, shapeLocations, shapeColors) {
      let r = 0.0;
      let g = 0.0;
      let b = 0.0;
      let a = 0.0;
      const x = this.thread.x;
      const y = this.thread.y;
      if (this.constants.pointsEnabled) {
        for (let i = 0; i < this.constants.pointCount; i++) {
          const px = pointsX[i];
          const py = pointsY[i];
          const dx = px - x;
          const dy = py - y;
          const d = Math.sqrt(dx * dx + dy * dy);
          const e = d / (this.constants.maxd * this.constants.pointCount);
          const power = this.constants.lightDecay;
          const e2 = Math.pow(this.constants.maxe, power) / (e == 0.0 ? 1.0 : Math.pow(e, power));
          r += pointsColorR[i] * e2;
          g += pointsColorG[i] * e2;
          b += pointsColorB[i] * e2;
        }
      }
      a = 1.0;
      const ws = this.constants.shapeWidth;
      const hs = this.constants.shapeHeight;
      for (let i = 0; i < this.constants.shapeCount; i++) {
        const posX = shapeLocations[i * 2];
        const posY = shapeLocations[i * 2 + 1];
        if (x >= posX && x < posX + ws && y >= posY && y < posY + hs) {
          const lx = x - posX;
          const ly = y - posY;
          const idx = 4 * (ly * ws + lx);
          const gidx = 4 * i * ws * hs + idx;
//          r += shapes[gidx + 0];
//          g += shapes[gidx + 1];
//          b += shapes[gidx + 2];
//          a += shapes[gidx + 3];
          const shapeSeed = i;
          const ar =
            createShapeGPU(
              this.constants.initRadius, this.constants.spikiness, shapeSeed, ws, hs, lx, ly);
          r += ar * shapeColors[i * 3 + 0];
          g += ar * shapeColors[i * 3 + 1];
          b += ar * shapeColors[i * 3 + 2];
          a += ar;
        }
      }
      return [r, g, b, 1];
    },
    {
      constants: {
        pointsEnabled: options.pointsEnabled,
        pointCount,
        maxd,
        maxe,
        lightDecay: options.lightDecay,
        shapeCount,
        shapeWidth: ws,
        shapeHeight: hs,
        initRadius: options.initRadius,
        spikiness: options.spikiness
      },
      output: [w, h]
    });

  const out =
    process(
      pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB,
      shapes, shapeLocations, shapeColors);

  // Get normalize factor
  let [maxr, maxg, maxb] = [0.0, 0.0, 0.0];
  let [avgr, avgg, avgb] = [0.0, 0.0, 0.0];
  const avg = (a, b) => a + b;
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const [r, g, b, a] = out[x][y];
      [maxr, maxg, maxb] = [Math.max(maxr, r), Math.max(maxg, g), Math.max(maxb, b)];
      [avgr, avgg, avgb] = [avg(avgr, r), avg(avgg, g), avg(avgb, b)];
    }
  }
  const m = w * h * 0.5;
  [avgr, avgg, avgb] = [m / (avgr == 0 ? m : avgr), m / (avgg == 0 ? m : avgg), m / (avgb == 0 ? m : avgb)];
  if (!options.pointsEnabled) {
    [avgr, avgg, avgb] = [1, 1, 1];
  }

  const f = (x) => Math.max(0, Math.min(255, Math.floor(x * 255)));
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const idx = 4 * (y * w + x);
      const [r, g, b, a] = out[x][y];
      data.data[idx + 0] = f(r * avgr);
      data.data[idx + 1] = f(g * avgg);
      data.data[idx + 2] = f(b * avgb);
      data.data[idx + 3] = f(a);
    }
  }

  ctx.putImageData(data, 0, 0);
};

const datgui = function() {
  const gui = new dat.GUI();
  gui.useLocalStorage = true;
  const defsize = Math.min(window.innerWidth, window.innerHeight) - 100;
  const options = {
    seed: 'hello',
    width: defsize,
    height: defsize,
    pointsEnabled: true,
    pointCount: 100,
    pointsRandomColor: true,
    pointsColor: [0, 0, 0, 0],
    lightDecay: 1.0,
    shapesEnabled: true,
    shapeWidth: 100,
    shapeHeight: 100,
    initRadius: 0.5,
    spikiness: 0.5
  };
  const reset = function() {
    main(options);
  };
  gui.remember(options);
  gui.add(options, 'seed').onFinishChange(reset);
  const size = gui.addFolder('Canvas Size');
  size.add(options, 'width', 0, 10000).onFinishChange(reset);
  size.add(options, 'height', 0, 10000).onFinishChange(reset);
  const points = gui.addFolder('Point Options');
  points.add(options, 'pointsEnabled').onFinishChange(reset);
  points.add(options, 'pointCount', 0, 1000).onFinishChange(reset);
  points.add(options, 'pointsRandomColor').onFinishChange(reset);
  points.addColor(options, 'pointsColor').onFinishChange(reset);
  const lights = gui.addFolder('Light Options');
  lights.add(options, 'lightDecay', 0.0, 10.0).onFinishChange(reset);
  const shapes = gui.addFolder('Shape Options');
  shapes.add(options, 'shapesEnabled').onFinishChange(reset);
  shapes.add(options, 'shapeWidth', 0.0, 500.0).onFinishChange(reset);
  shapes.add(options, 'shapeHeight', 0.0, 500.0).onFinishChange(reset);
  shapes.add(options, 'initRadius', 0.0, 5.0).onFinishChange(reset);
  shapes.add(options, 'spikiness', 0.0, 5.0).onFinishChange(reset);
  reset();
};

//$(main({width: 400, height: 400}));
$(datgui);

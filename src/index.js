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

  GPU.addNativeFunction('rand',
  `
	float rand(vec2 co){
			return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
	}
	`);

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
			.css('position', 'absolute')
      .css('top', '0px')
      .css('z-index', '-1')
      .css('width', '100%')
      .css('height', '100vh')
      .css('background', 'black')
      .css('display', 'flex')
      .css('justify-content', 'center')
      .css('align-items', 'center')
      .append(canvas);
  }

  $("#mainImage")
    .css('border', '1px solid white')

	const ws = options.shapeWidth;
	const hs = options.shapeHeight;
  const shapeCount = options.shapeCount;
  const shapeVariations = options.shapeVariations;
  const shapeLocations = new Float32Array(shapeCount * 2);
  const shapeColors = new Float32Array(shapeCount * 3);
  const shapeSizes = new Float32Array(shapeCount);
  const maxShapesPerPixel = 10;
  const shapeAllocations = new Float32Array(Math.ceil(maxShapesPerPixel * Math.ceil(w / ws) * Math.ceil(h / hs)));
  for (let i = 0; i < shapeCount; i++) {
    const color = {
      r: options.shapesRandomColor ? Math.random() : options.shapesColor[0],
      g: options.shapesRandomColor ? Math.random() : options.shapesColor[1],
      b: options.shapesRandomColor ? Math.random() : options.shapesColor[2],
    };
    shapeColors[i * 3 + 0] = color.r;
    shapeColors[i * 3 + 1] = color.g;
    shapeColors[i * 3 + 2] = color.b;
    shapeSizes[i] = 1.0;
    if (options.shapeDistribution === 'Circles') {
      const circleCount = options.shapeCircleCount;
      const perCircle = Math.floor(shapeCount / circleCount);
      const circle = Math.floor(i / perCircle) + 1;
      const r = options.shapeCircleRadius * circle;
      const angle = 2.0 * Math.PI * (i - circle * perCircle) / perCircle;
      const x = Math.floor(0.5 * w + r * Math.cos(angle));
      const y = Math.floor(0.5 * h + r * Math.sin(angle));
      shapeLocations[i * 2] = x;
      shapeLocations[i * 2 + 1] = y;
      shapeSizes[i] /= (circleCount - circle + 2);
    } else if (options.shapeDistribution === 'UniformRandomOffset') {
      const sq = Math.ceil(Math.sqrt(shapeCount));
      const sq2 = Math.floor(Math.sqrt(shapeCount));
      const wstep = Math.floor(w / sq2);
      const hstep = Math.floor(h / sq2);
      const rstep = options.shapeCircleRadius;
      const x = Math.floor(i / sq) * hstep + rstep * (Math.random() - 0.5) + 0.5 * ws;
      const y = Math.floor(i % sq) * wstep + rstep * (Math.random() - 0.5) + 0.5 * hs;
      shapeLocations[i * 2] = x;
      shapeLocations[i * 2 + 1] = y;
    } else {
      shapeLocations[i * 2] = Math.floor(Math.random() * w);
      shapeLocations[i * 2 + 1] = Math.floor(Math.random() * h);
    }
    const xr = Math.floor(shapeLocations[i * 2] / ws);
    const wr = ws * shapeSizes[i];
    const yr = Math.floor(shapeLocations[i * 2 + 1] / hs);
    const hr = hs * shapeSizes[i];
    const wr5 = Math.floor(0.5 * wr);
    const hr5 = Math.floor(0.5 * hr);
    {
      const k = 1;
      for (let x = -k; x <= k; x += 1) {
        for (let y = -k; y <= k; y += 1) {
          let idx = maxShapesPerPixel * ((yr + y) * Math.ceil(w / ws) + (xr + x));
          let shapeIndex = 0;
          while (shapeAllocations[idx] > 0 && shapeIndex < maxShapesPerPixel) {
            idx += maxShapesPerPixel;
            shapeIndex++;
          }
          shapeAllocations[idx] = i + 1;
        }
      }
    }
  }

  const gpu = new GPU.GPU({ canvas: canvas });
	addNoiseFunctions(gpu);
  function lerp2(x, y, a) {
    return x * (1.0 - a) + y * a;
  }
  function snoise3(a, b, c) {
    const ret = snoise([a, b, c]);
    return ret;
  }
  function createShapeGPU(spikeFrequency, initRadius, spikiness, shapeSeed, ws, hs, x, y) {
    const dx = x;
    const dy = y;
    const radius = 0.5 * Math.max(ws, hs);
    const vradius =
      (0.75 + 0.25 * Math.sin(10.0 * Math.atan2(dy, dx))) * radius;
    const angle = Math.atan2(dy, dx) + Math.PI;
    const nangle = angle / Math.PI / 2.0;
    const c = spikeFrequency;
    const rnoise = 0.5 * (snoise3(Math.sin(c * angle), shapeSeed, 0.0) + 1.0);
    const vradius2 = (initRadius + spikiness * rnoise) * radius;
    const currentRadius = Math.sqrt(dx * dx + dy * dy);
    const border = currentRadius >= vradius2 && currentRadius <= vradius2 + 1.0;
    const bc1 = snoise3((0.5 * x) / ws, (0.5 * y) / ws, shapeSeed);
    const bc2 = 0.5 * bc1 + 1.0;
          
    // Light emission from shape
    const da = Math.max(1.0, Math.abs(currentRadius - vradius2));
    const d = currentRadius <= vradius2 && this.constants.shapeInteriors ? 1.0 : da;
    const e = d / (this.constants.maxd * this.constants.shapeCount);
    const power = this.constants.lightDecay;
    const light = 1.0 * Math.pow(this.constants.maxe, power) / (e == 0.0 ? 1.0 : Math.pow(e, power));
    const lightedShape = light;

    const shape = currentRadius <= vradius2 ? lightedShape : 1.0 * lightedShape;

    return shape;
  };
	gpu.setFunctions([lerp2, snoise3, createShapeGPU]);

  const maxd = Math.sqrt(w * w + h * h);
  const maxe = 1.0 / shapeCount;
  const process = gpu
    .createKernel(function (
      iteration, out, //pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB,
      shapeLocations, shapeColors, shapeSizes, shapeAllocations) {
      const x = this.thread.x;
      const y = this.thread.y;
      const prev = out[y][x];
      let r = prev[0];
      let g = prev[1];
      let b = prev[2];
      let a = prev[3];
//      if (this.constants.pointsEnabled) {
//        for (let i = 0; i < this.constants.pointCount; i++) {
//          const px = pointsX[i];
//          const py = pointsY[i];
//          const dx = px - x;
//          const dy = py - y;
//          const d = Math.sqrt(dx * dx + dy * dy);
//          const e = d / (this.constants.maxd * this.constants.pointCount);
//          const power = this.constants.lightDecay;
//          const e2 = Math.pow(this.constants.maxe, power) / (e == 0.0 ? 1.0 : Math.pow(e, power));
//          r += pointsColorR[i] * e2;
//          g += pointsColorG[i] * e2;
//          b += pointsColorB[i] * e2;
//        }
//      }
      a = 1.0;
      if (this.constants.shapesEnabled) {
        const ws = this.constants.shapeWidth;
        const hs = this.constants.shapeHeight;
        const i = iteration;
        const shapeVariationIndex = i % this.constants.shapeVariations;
        const posX = shapeLocations[i * 2];
        const posY = shapeLocations[i * 2 + 1];
        const size = shapeSizes[i];
        const halfws = Math.floor(0.5 * ws * size);
        const halfhs = Math.floor(0.5 * hs * size);
        const lx = (x - posX) * ws / hs;
        const ly = (y - posY) * hs / ws;
        const shapeSeed = shapeVariationIndex;
        const initRadiusMul = this.constants.randomizeInitRadius ? rand([5 * i, 19 * i]) : 1.0;
        const spikinessMul = this.constants.randomizeSpikiness ? rand([13 * i, 7 * i]) : 1.0;
        const ar =
          createShapeGPU(
            this.constants.spikeFrequency,
            initRadiusMul * this.constants.initRadius,
            spikinessMul * this.constants.spikiness,
            shapeSeed,
            ws * size,
            hs * size,
            lx,
            ly);
        r += ar * shapeColors[i * 3 + 0];
        g += ar * shapeColors[i * 3 + 1];
        b += ar * shapeColors[i * 3 + 2];
        a += ar;
      }
      return [r, g, b, 1];
    },
    {
      constants: {
        width: w,
        height: h,
        shapesEnabled: options.shapesEnabled,
        shapeCount: options.shapeCount,
        maxd,
        maxe,
        lightDecay: options.lightDecay,
        shapeCount: 1,
        shapeVariations,
        shapeWidth: ws,
        shapeHeight: hs,
        initRadius: options.initRadius,
        spikiness: options.spikiness,
        randomizeInitRadius: options.randomizeInitRadius,
        randomizeSpikiness: options.randomizeSpikiness,
        spikeFrequency: options.spikeFrequency,
        shapeInteriors: options.shapeInteriors,
        maxShapesPerPixel,
      },
      output: [w, h],
      immutable: true,
      pipeline: true
    }).setPipeline(true);

  const before = new Date().getTime();
  let temp = gpu.createKernel(function() { return [0, 0, 0, 0]; }).setPipeline(true).setOutput([w, h])();
  for (let i = 0; i < shapeCount; i++) {
    const temp2 =
      process(i, temp,
        //pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB,
        shapeLocations, shapeColors, shapeSizes, shapeAllocations);
    temp.delete();
    temp = temp2;
  }
	
	const maxrows = gpu.createKernel(function(image) {
			let pixel = [0, 0, 0, 0];
			for (let x = 0; x < this.constants.width; x++) {
				const p2 = image[this.thread.y][x];
		    const r = Math.max(pixel[0], p2[0]);	
		    const g = Math.max(pixel[1], p2[1]);
		    const b = Math.max(pixel[2], p2[2]);	
		    const a = Math.max(pixel[3], p2[3]);	
        pixel = [r, g, b, a];
			}
			return pixel;
		})
    .setConstants({
      width: w,
      height: h
    })
    .setPipeline(true)
 		.setOutput([h]);
	const maxrowsOutput = maxrows(temp);

  const maxcols = gpu.createKernel(function(image) {
			let pixel = [0, 0, 0, 0];
			for (let y = 0; y < this.constants.height; y++) {
				const p2 = image[y];
		    const r = Math.max(pixel[0], p2[0]);	
		    const g = Math.max(pixel[1], p2[1]);
		    const b = Math.max(pixel[2], p2[2]);	
		    const a = Math.max(pixel[3], p2[3]);	
        pixel = [r, g, b, a];
			}
			const m1 = Math.max(pixel[0], pixel[1]);
      const m2 = Math.max(pixel[2], m1);
      return 1.0e-1 / m2;
		})
    .setConstants({
      width: w,
      height: h
    })
    .setPipeline(true)
 		.setOutput([1]);
	const normalize = maxcols(maxrowsOutput);

	const kernel = gpu.createKernel(function(image, brightness, normalize) {
			const pixel = image[this.thread.y][this.thread.x];
      const n = brightness * normalize[0];
			this.color(pixel[0] * n, pixel[1] * n, pixel[2] * n, pixel[3] * n);
		})
		.setGraphical(true)
 		.setOutput([w, h]);
	kernel(temp, options.brightness, normalize);
  //ctx.drawImage(kernel.getCanvas(), 0, 0)
  // $("body").append(kernel.canvas);

  const after = new Date().getTime();
  const elapsed = after - before;
  console.log("Time elapsed:", elapsed / 1000.0);

  //const out = temp.toArray();
  temp.delete();

  // Get normalize factor
//  let [maxr, maxg, maxb] = [0.0, 0.0, 0.0];
//  let [avgr, avgg, avgb] = [0.0, 0.0, 0.0];
//  const avg = (a, b) => a + b;
//  const rs = [];
//  const gs = [];
//  const bs = [];
//  for (let x = 0; x < w; x++) {
//    for (let y = 0; y < h; y++) {
//      const [r, g, b, a] = out[x][y];
//      [maxr, maxg, maxb] = [Math.max(maxr, r), Math.max(maxg, g), Math.max(maxb, b)];
//      [avgr, avgg, avgb] = [avg(avgr, r), avg(avgg, g), avg(avgb, b)];
//      rs.push(r);
//      gs.push(g);
//      bs.push(b);
//    }
//  }
//  rs.sort();
//  gs.sort();
//  bs.sort();
//  const meanr = rs[Math.floor(rs.length * 0.5)];
//  const meang = gs[Math.floor(gs.length * 0.5)];
//  const meanb = bs[Math.floor(bs.length * 0.5)];

//  let [stddevr, stddevg, stddevb] = [0.0, 0.0, 0.0];
//  const sq = x => x * x;
//  const lbound = Math.floor(rs.length * 0.0);
//  const ubound = Math.floor(rs.length * 1.0 - 1.0e-3);
//  for (let i = lbound; i < ubound; i++) {
//    stddevr += sq(rs[i] - meanr);
//    stddevg += sq(gs[i] - meang);
//    stddevb += sq(bs[i] - meanb);
//  }
//  stddevr = Math.sqrt(stddevr);
//  stddevg = Math.sqrt(stddevg);
//  stddevb = Math.sqrt(stddevb);

//  const lr = 0;
//  const lg = 0;
//  const lb = 0;
//  const mm = Math.max(Math.max(rs[ubound], gs[ubound]), bs[ubound]);
//  const mr = mm;
//  const mg = mm;
//  const mb = mm;
//  const dr = mr - lr;
//  const dg = mg - lg;
//  const db = mb - lb;

//  const f = (x) => Math.max(0, Math.min(255, Math.floor(options.brightness * x * 255)));
//  const ex = x => x;
//  for (let x = 0; x < w; x++) {
//    for (let y = 0; y < h; y++) {
//      const idx = 4 * (y * w + x);
//      const [r, g, b, a] = out[x][y];
//      const r2 = (r - lr) / dr;
//      const g2 = (g - lg) / dg;
//      const b2 = (b - lb) / db;
//      data.data[idx + 0] = f(r2);
//      data.data[idx + 1] = f(g2);
//      data.data[idx + 2] = f(b2);
//      data.data[idx + 3] = f(a);
//    }
//  }

  //ctx.putImageData(data, 0, 0);
};

const datgui = function() {
  const gui = new dat.GUI();
  gui.useLocalStorage = true;
  const defsize = Math.min(window.innerWidth, window.innerHeight) - 100;
  const options = {
    seed: 'hello2',
    width: defsize,
    height: defsize,
    brightness: 1,
    lightDecay: 1.2,
    shapesEnabled: true,
    shapeInteriors: false,
    shapesRandomColor: true,
    shapesColor: [0, 0, 0, 0],
    shapeCount: 100,
    shapeVariations: 10,
    shapeDistribution: 'Circular',
    shapeCircleCount: 5,
    shapeCircleRadius: 50.0,
    shapeWidth: 200,
    shapeHeight: 200,
    initRadius: 0.5,
    randomizeInitRadius: false,
    spikiness: 0.5,
    randomizeSpikiness: false,
		spikeFrequency: 8.0
  };
  const reset = function() {
    main(options);
  };
  gui.width = 300;
  gui.remember(options);
  gui.add(options, 'seed').onFinishChange(reset);
  const size = gui.addFolder('Canvas');
  size.add(options, 'width', 0, 10000).onFinishChange(reset);
  size.add(options, 'height', 0, 10000).onFinishChange(reset);
  size.add(options, 'brightness', 0, 10000).onFinishChange(reset);
  const lights = gui.addFolder('Light Options');
  lights.add(options, 'lightDecay', 0.0, 10.0).onFinishChange(reset);
  const shapes = gui.addFolder('Shape Options');
  shapes.add(options, 'shapesEnabled').onFinishChange(reset);
  shapes.add(options, 'shapeInteriors').onFinishChange(reset);
  shapes.add(options, 'shapesRandomColor').onFinishChange(reset);
  shapes.addColor(options, 'shapesColor').onFinishChange(reset);
  shapes.add(options, 'shapeWidth', 0.0, 500.0).onFinishChange(reset);
  shapes.add(options, 'shapeHeight', 0.0, 500.0).onFinishChange(reset);
  shapes.add(options, 'shapeCount', 0, 1000).onFinishChange(reset);
  shapes.add(options, 'shapeVariations', 0, 1000).onFinishChange(reset);
  shapes.add(options, 'shapeDistribution',
    ['Random', 'UniformRandomOffset', 'Circles']).onFinishChange(reset);
  shapes.add(options, 'shapeCircleCount', 0, 20).onFinishChange(reset);
  shapes.add(options, 'shapeCircleRadius', 0, 500).onFinishChange(reset);
  shapes.add(options, 'initRadius', 0.0, 5.0).onFinishChange(reset);
  shapes.add(options, 'randomizeInitRadius').onFinishChange(reset);
  shapes.add(options, 'spikiness', 0.0, 5.0).onFinishChange(reset);
  shapes.add(options, 'randomizeSpikiness').onFinishChange(reset);
  shapes.add(options, 'spikeFrequency', 0.0, 10.0).onFinishChange(reset);
  reset();
};

//$(main({width: 400, height: 400}));
$(datgui);

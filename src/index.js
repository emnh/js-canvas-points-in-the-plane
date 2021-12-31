const $ = require('jquery');
const dat = require('dat.gui');
import GPU from 'gpu.js';
import SimplexNoise from 'simplex-noise';

const seedrandom = require('seedrandom');

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
      const innerCircle = Math.sqrt(dx * dx + dy * dy) <= radius;
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

	const ws = 100;
	const hs = 100;
  const shapeCount = 100;	
  const shapeLocations = new Float32Array(shapeCount * 2);
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
    for (let x = 0; x < ws; x++) {
      for (let y = 0; y < hs; y++) {
        const idx = 4 * (y * ws + x);
        const gidx = i * ws * hs + idx;
        shapes[gidx + 0] = shape[idx + 0] * color.r;
        shapes[gidx + 1] = shape[idx + 1] * color.g;
        shapes[gidx + 2] = shape[idx + 2] * color.b;
        shapes[gidx + 3] = shape[idx + 3];
      }
    }
	}


  const gpu = new GPU.GPU();
  const maxd = Math.sqrt(w * w + h * h);
  const maxe = 1.0 / pointCount;
  const process = gpu
    .createKernel(function (pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB, shapes, shapeLocations) {
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
          const idx = 4 * ((y - posY) * ws + (x - posX));
          const gidx = i * ws * hs + idx;
          r += shapes[gidx + 0];
          g += shapes[gidx + 1];
          b += shapes[gidx + 2];
          a += shapes[gidx + 3];
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
        shapeHeight: hs
      },
      output: [w, h]
    });

  const out = process(pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB, shapes, shapeLocations);

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
    lightDecay: 1.0
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
  reset();
};

//$(main({width: 400, height: 400}));
$(datgui);

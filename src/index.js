const $ = require('jquery');
const dat = require('dat.gui');
import GPU from 'gpu.js';

const main = function(options) {

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
  const pointCount = 100;
  const pointsX = new Float32Array(pointCount);
  const pointsY = new Float32Array(pointCount);
  const pointsColorR = new Float32Array(pointCount);
  const pointsColorG = new Float32Array(pointCount);
  const pointsColorB = new Float32Array(pointCount);


  for (let i = 0; i < pointCount; i++) {
    const c = {
      r: Math.random(),
      g: Math.random(),
      b: Math.random(),
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


  const gpu = new GPU.GPU();
  const maxd = Math.sqrt(w * w + h * h);
  const maxe = 1.0 / pointCount;
  const process = gpu
    .createKernel(function (pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB) {
      let r = 0.0;
      let g = 0.0;
      let b = 0.0;
      const x = this.thread.x;
      const y = this.thread.y;
      for (let i = 0; i < this.constants.pointCount; i++) {
        const px = pointsX[i];
        const py = pointsY[i];
        const dx = px - x;
        const dy = py - y;
        const d = Math.sqrt(dx * dx + dy * dy);
        const e = d / (this.constants.maxd * this.constants.pointCount);
        const power = 1.0;
        const e2 = Math.pow(this.constants.maxe, power) / (e == 0.0 ? 1.0 : Math.pow(e, power));
        r += pointsColorR[i] * e2;
        g += pointsColorG[i] * e2;
        b += pointsColorB[i] * e2;
      }
      return [r, g, b, 1];
    },
    {
      constants: {
        pointCount,
        maxd,
        maxe,
      },
      output: [w, h]
    });

  const out = process(pointsX, pointsY, pointsColorR, pointsColorG, pointsColorB);

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
    width: defsize,
    height: defsize
  };
  const reset = function() {
    main(options);
  };
  const size = gui.addFolder('Canvas Size');
  size.add(options, 'width', 0, 10000).onFinishChange(reset);
  size.add(options, 'height', 0, 10000).onFinishChange(reset);
  reset();
};

//$(main({width: 400, height: 400}));
$(datgui);

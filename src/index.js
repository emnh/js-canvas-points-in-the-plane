const $ = require('jquery');
import GPU from 'gpu.js';

const main = function() {
  const w = 800;
  const h = 800;

  const canvas = document.createElement('canvas');

  canvas.width = w;
  canvas.height = h;

  document.body.appendChild(canvas);

  const ctx = canvas.getContext('2d');

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
      x: Math.random() * w,
      y: Math.random() * h,
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
        const e2 = Math.pow(this.constants.maxe, power) / Math.pow(e, power);
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
  [avgr, avgg, avgb] = [m / avgr, m / avgg, m / avgb];

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

$(main);

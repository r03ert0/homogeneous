"use strict";
const cli = require('cli');
const PLY = require('ply.js');
var iter = 0;
cli.parse({
    input: ['i', 'Path to mesh to process, in PLY format', 'file'],
    niter: ['n', 'Number of iterations, default 5000', 'int', 5000],
    alpha: ['s', 'Displacement fraction, from 0 to 0.5, default 0.3', 'float', 0.3],
    tau: ['t', 'Displacement smoothness, from 0 to 1, default 0.1', 'float', 0.1],
    aspectRatio: ['a', 'Maximum aspect ratio, 1 is an equilateral triangle, the default is 25 which is quite ugly', 'float', 25],
    output: ['o', 'Path to directory resulting mesh, in PLY format', 'file', 'output.ply']
});
var input, output, niter, alpha, tau, aspectRatio;
cli.main(function (args, options) {
    if (options.input) {
        console.log('input',options.input);
        input = options.input;
    }
    if (options.output) {
        console.log('output',options.output);
        output = options.output;
    }
    if (options.niter) {
        console.log('niter',options.niter);
        niter = options.niter;
    }
    if (options.alpha) {
        console.log('alpha',options.alpha);
        alpha = options.alpha;
    }
    if (options.tau) {
        console.log('tau',options.tau);
        tau = options.tau;
    }
    if (options.aspectRatio) {
        console.log('aspectRatio',options.aspectRatio);
        aspectRatio = options.aspectRatio;
    }
});

var params = new (function() {
    this.alpha = alpha;
    this.tau = tau;
    this.maxAspectRatio = aspectRatio;
})();


//const name = 'icosahedron.ply';
//const name = 'surf.sphere-502.ply';
//const name = 'surf.sphere-2502.ply';
//const name = 'surf.sphere-14855-laplace5k.ply';

const mesh = PLY.loadPLY(input);
for(iter=0;iter<niter;iter++) {
    step();
}
PLY.savePLY(mesh, output);
console.log('done');

function add3D(a, b) {
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]];
}
function sub3D(a, b) {
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]];
}
function sca3D(a, t) {
    return [a[0]*t,a[1]*t,a[2]*t];
}
function dot3D(a, b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
function cross3D(a, b) {
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
}

function norm3D(v) {
    return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
function triangle_area(p0, p1, p2)
{
    let a,b,c; // side lengths
    let s; // semiperimeter
    let area;

    a=norm3D(sub3D(p0,p1));
    b=norm3D(sub3D(p1,p2));
    c=norm3D(sub3D(p2,p0));
    s=(a+b+c)/2.0;

    if(s*(s-a)*(s-b)*(s-c)<0)
        area=0;
    else
        area=Math.sqrt(s*(s-a)*(s-b)*(s-c));
    
    return area;
}
function step() {
    if(params.stopSimulation) {
        return;
    }

    const [p, t] = [mesh.p, mesh.t];
    const np = p.length;
    const nt = t.length;
    const f = [];
    const n = [];
    const deg = [];
    const min = [];
    let A;
    let ia, ib, ic;
    let a, b, c;
    let ab, bc, ca;
    let v;
    let dab, dbc, dca;
    let ha, hb, hc;
    let i, num;
    const beta = [];
    let x;
    let s, ar;
    let fa, fb, fc;

    // initialise force
    for(i=0;i<np;i++) {
        f[i] = [0,0,0];
        n[i] = 0;
        beta[i] = Infinity;
    }

    for(i=0;i<nt;i++) {
        deg[i] = 0;
    }

    // compute repulsion forces and find minimum triangle distance
    for(i=0;i<nt;i++) {
        // get vertex indices
        [ia, ib, ic]=[...t[i]];

        // get vertices
        a=[...p[ia]];
        b=[...p[ib]];
        c=[...p[ic]];

        // triangle area
        A=triangle_area(a,b,c);

        // edge lengths
        [ab, bc, ca]=[sub3D(a, b), sub3D(b, c), sub3D(c, a)];
        [dab, dbc, dca] = [norm3D(ab), norm3D(bc), norm3D(ca)];

        // triangle aspect ratio
        s = (dab+dbc+dca)/2;
        ar = dab*dbc*dca/(8*(s-dab)*(s-dbc)*(s-dca));

        // mark degenerate triangles
        if(ar>params.maxAspectRatio) {
            deg[i] = ar;
        }

        // triangle altitudes
        [ha, hb, hc] = [A/dbc, A/dca, A/dab];

        // compute the minimum triangle dimension
        min[i] = Math.min(dab, dbc, dca, ha, hb, hc);

        // repulsion forces
        // F = bcx(abxbc)xcb/|bc|^2
        fa = sca3D(cross3D(bc,cross3D(ab,bc)),1/dbc/dbc);
        fb = sca3D(cross3D(ca,cross3D(bc,ca)),1/dca/dca);
        fc = sca3D(cross3D(ab,cross3D(ca,ab)),1/dab/dab);
        fa = sca3D(fa,1/Math.pow(ha,3));
        fb = sca3D(fb,1/Math.pow(hb,3));
        fc = sca3D(fc,1/Math.pow(hc,3));

        // add forces
        f[ia]=add3D(f[ia], fa);
        f[ib]=add3D(f[ib], fb);
        f[ic]=add3D(f[ic], fc);
        n[ia]++;
        n[ib]++;
        n[ic]++;
    }
    for(i=0;i<np;i++) {
        f[i]=sca3D(f[i],1/n[i]);
    }

    let max, tmp;
    for(i=0;i<nt;i++) {
        if(deg[i] > 0) {
            continue;
        }
        max = Math.max(norm3D(f[t[i][0]]),norm3D(f[t[i][1]]),norm3D(f[t[i][2]]));
        beta[t[i][0]] = Math.min(beta[t[i][0]], params.alpha*min[i]/max);
        beta[t[i][1]] = Math.min(beta[t[i][1]], params.alpha*min[i]/max);
        beta[t[i][2]] = Math.min(beta[t[i][2]], params.alpha*min[i]/max);
    }
    for(i=0;i<np;i++) {
        if(beta[i]>1) {
            beta[i]=1;
        }
    }
    for(i=0;i<np;i++) {
        f[i]=sca3D(f[i], beta[i]);
    }

    // prevent triangles to flip, and degenerate triangles to worsen
    let at, bt, ct;
    num = 0;
    for(i=0;i<nt;i++) {
        [ia, ib, ic]=[...t[i]];
        a=p[ia];
        b=p[ib];
        c=p[ic];
        at=add3D(a, f[ia]);
        bt=add3D(b, f[ib]);
        ct=add3D(c, f[ic]);
        at = sca3D(at,1/norm3D(at));
        bt = sca3D(bt,1/norm3D(bt));
        ct = sca3D(ct,1/norm3D(ct));
        ab=sub3D(at, bt);
        bc=sub3D(bt, ct);
        ca=sub3D(ct, at);
        if(deg[i] > 0) {
            dab=norm3D(ab);
            dbc=norm3D(bc);
            dca=norm3D(ca);
            s = (dab+dbc+dca)/2;
            ar = dab*dbc*dca/(8*(s-dab)*(s-dbc)*(s-dca));
            // worsening degenerate triangle
            if(ar > deg[i]) {
                f[ia]=[0,0,0];
                f[ib]=[0,0,0];
                f[ic]=[0,0,0];
                num++;
            }
        }
        // inverting triangle
        if(dot3D(cross3D(ab,bc),a) < 0) {
            f[ia]=[0,0,0];
            f[ib]=[0,0,0];
            f[ic]=[0,0,0];
            num++;
        }
    }

    // apply displacement to mesh
    let totDisp = 0;
    let g = params.tau
    for(i=0;i<np;i++) {
        v=add3D(p[i],f[i]);
        v=sca3D(v,1/norm3D(v));
        p[i] = add3D(sca3D(p[i],1-g), sca3D(v,g));
        totDisp += norm3D(sub3D(v, p[i]));
    }

    if(iter%100 === 0) {
        console.log(`${iter} total displacement: ${totDisp}`);
    }
}

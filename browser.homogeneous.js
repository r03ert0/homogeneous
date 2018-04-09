"use strict";
var scene,renderer,camera,cameraControls, mesh, forceDisplay;
var dataset, labels=[], Y
var iter = 0;
var params = new (function() {
    this.message = 'Homogeneous';
    this.alpha = 0.3;
    this.tau = 0.1;
    this.maxAspectRatio = 25;
    this.showForces = true;
    this.showVertexIndices = false;
    this.stopSimulation = false;
    this.savePLY= savePLY;
})();
var gui = new dat.GUI();
gui.add(params, 'message');
gui.add(params, 'alpha');
gui.add(params, 'tau');
gui.add(params, 'maxAspectRatio');
gui.add(params, 'showForces');
gui.add(params, 'showVertexIndices');
gui.add(params, 'stopSimulation');
gui.add(params, 'savePLY');

init();
loadMesh();

function init() {
    var w, h;
    var container=$('#view');
    w = container.width();
    h = container.height();

// init the scene
    renderer = new THREE.WebGLRenderer({
        antialias: true,    // to get smoother output
        preserveDrawingBuffer: true    // to allow screenshot
    });
    renderer.setSize( w, h );
    renderer.setPixelRatio(window.devicePixelRatio);
    container[0].appendChild(renderer.domElement);
    scene = new THREE.Scene();
    scene.background = new THREE.Color( 0xffffff );
    camera    = new THREE.PerspectiveCamera(50, w/h, 1, 100 );
    camera.position.set(0, 0, 5);
    scene.add(camera);

    // create a camera control
    cameraControls=new THREE.TrackballControls(camera,container[0] )
    cameraControls.rotateSpeed = 6;
}

function initForceDisplay() {
    const p = mesh.geometry.attributes.position.array;
    const np = p.length/3;
    const geometry=new THREE.Geometry();
    let i;

    for(i=0;i<np;i++) {
        geometry.vertices.push(new THREE.Vector3(0,0,0));
        geometry.vertices.push(new THREE.Vector3(0,0,0));
    }
    var material=new THREE.LineBasicMaterial({color:'black'});
    forceDisplay=new THREE.LineSegments(geometry,material);
    scene.add(forceDisplay);
}

function loadMesh() {
    var loader = new THREE.PLYLoader();
//    const name = 'icosahedron.ply';
//    const name = 'surf.sphere-502.ply';
//    const name = 'surf.sphere-2502.ply';
    const name = 'surf.sphere-14855-laplace5k.ply';
    loader.load( name, function ( geometry ) {
        geometry.computeVertexNormals();
        var material = new THREE.MeshNormalMaterial({wireframe: true});
        mesh = new THREE.Mesh( geometry, material );
        mesh.castShadow = true;
        mesh.receiveShadow = true;
        scene.add( mesh );

        // deform
        //deform();

        // init force display
        initForceDisplay();

        // trigger animation
        animate();
    } );
}

function deform() {
    const c = [0,0,0];
    let i,j,np,v;
    const p = mesh.geometry.attributes.position.array;
    np = p.length/3;
    for(i=0;i<np;i++) {
        v = [p[3*i+0]-c[0],p[3*i+1]-c[1],p[3*i+2]-c[2]];
        length = Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        p[3*i+0]=v[0]/length;
        p[3*i+1]=v[1]/length;
        p[3*i+2]=v[2]/length;

        $('body').append(`<div class='vertex' id='v${i}'>${i}</div>`);
    }
}

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

    const p = mesh.geometry.attributes.position.array;
    const t = mesh.geometry.index.array;
    const np = p.length/3;
    const nt = t.length/3;
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
        [ia, ib, ic]=[t[3*i+0], t[3*i+1], t[3*i+2]];

        // get vertices
        a=[p[3*ia+0], p[3*ia+1], p[3*ia+2]];
        b=[p[3*ib+0], p[3*ib+1], p[3*ib+2]];
        c=[p[3*ic+0], p[3*ic+1], p[3*ic+2]];

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
        max = Math.max(norm3D(f[t[3*i+0]]),norm3D(f[t[3*i+1]]),norm3D(f[t[3*i+2]]));
        beta[t[3*i+0]] = Math.min(beta[t[3*i+0]], params.alpha*min[i]/max);
        beta[t[3*i+1]] = Math.min(beta[t[3*i+1]], params.alpha*min[i]/max);
        beta[t[3*i+2]] = Math.min(beta[t[3*i+2]], params.alpha*min[i]/max);
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
        [ia, ib, ic]=[t[3*i+0], t[3*i+1], t[3*i+2]];
        a=[p[3*ia+0], p[3*ia+1], p[3*ia+2]];
        b=[p[3*ib+0], p[3*ib+1], p[3*ib+2]];
        c=[p[3*ic+0], p[3*ic+1], p[3*ic+2]];
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
        v=add3D([p[3*i+0],p[3*i+1],p[3*i+2]],f[i]);
        v=sca3D(v,1/norm3D(v));
        [p[3*i+0],p[3*i+1],p[3*i+2]] = [p[3*i+0]*(1-g)+v[0]*g,p[3*i+1]*(1-g)+v[1]*g,p[3*i+2]*(1-g)+v[2]*g];
        totDisp += norm3D(sub3D(v, [p[3*i+0],p[3*i+1],p[3*i+2]]));
    }

    // draw forces
    if(params.showForces) {
        for(i=0;i<np;i++) {
            v = [p[3*i+0],p[3*i+1],p[3*i+2]];
            forceDisplay.geometry.vertices[2*i].x=v[0];
            forceDisplay.geometry.vertices[2*i].y=v[1];
            forceDisplay.geometry.vertices[2*i].z=v[2];
            forceDisplay.geometry.vertices[2*i+1].x=v[0]+f[i][0]*100;
            forceDisplay.geometry.vertices[2*i+1].y=v[1]+f[i][1]*100;
            forceDisplay.geometry.vertices[2*i+1].z=v[2]+f[i][2]*100;
        }
        forceDisplay.geometry.verticesNeedUpdate=true;
    }

    if(iter%100 === 0) {
        console.log(`${iter} total displacement: ${totDisp}`);
    }
}

function animate() {
    step();
    mesh.geometry.attributes.position.needsUpdate = true;
    iter++;

    if(typeof mesh !== 'undefined' && params.showVertexIndices) {
        var i,p,P=mesh.geometry.attributes.position.array,F=mesh.geometry.index.array;
        // show vertex indices
        for(i=0;i<P.length/3;i++) {
            p=screenXY(P[3*i+0],P[3*i+1],P[3*i+2]);
            $("#v"+i).css({left:p.x,top:p.y});
        }
    }

    requestAnimationFrame( animate );
    render();
}
function render() {
    cameraControls.update();
    renderer.render( scene, camera );
}
/**
 * @function screenXY
 */
function screenXY(x,y,z) {
    var vector = new THREE.Vector3(x,y,z);
    var canvas = $('canvas')[0];

    // map to normalized device coordinate (NDC) space
    vector.project( camera );

    // map to 2D screen space
    return {
        x: Math.round( (   vector.x + 0.9 ) * canvas.width/ 4 ),
        y: Math.round( ( - vector.y + 0.9 ) * canvas.height / 4 )
    }
}

function savePLY() {
    const parr = mesh.geometry.attributes.position.array;
    const tarr = mesh.geometry.index.array;
    const np = parr.length/3;
    const nt = tarr.length/3;
    const p = [];
    const t = [];
    let i;
    for(i=0;i<np;i++) {
        p.push(parr.slice(3*i,3*i+3));
    }
    for(i=0;i<nt;i++) {
        t.push(tarr.slice(3*i,3*i+3));
    }
    const str = PLY.encodePLY({p:p, t:t});

    const a = document.createElement('a');
    a.href = 'data:text/csv;charset=utf-8,' + str;
    let name = prompt("Save PLY mesh As...", "mesh.ply");
    if(name !== null) {
        a.download=name;
        document.body.appendChild(a);
        a.click();
    }
}

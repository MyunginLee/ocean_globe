// earth . 
#include "headers_10.hpp"

using namespace al;
using namespace std;
using namespace gam;

Mesh cell_sphere, earth_sphere, back_mesh, spike_mesh;
Pan<> mPan;
Env<3> mAmpEnv;
float OSCcarMul, OSCmodMul, OSCvibRate, OSCvibDepth, OSCtable,tableH;

LFO<> mod;
OnePole<> onePole;
string slurp(string fileName);
Vec3f randomVec3f(float scale);
static bool fullscreen = false;

struct Cells
{
  Vec3f pos, vel, acc, axis, vib, ang_moment;
  Color col;
  float angle;
  void update(float dt)
  {
    pos += vel * dt;
    // float dt = dt_ms;
  }
  void draw(Graphics &g, Mesh &m)
  {
    g.pushMatrix();
    g.translate(pos);
    g.rotate(angle, vel);
    g.color(col);
    g.draw(cell_sphere);
    g.popMatrix();
  }
};

struct earths
{
  Vec3f pos, vel, acc, axis, vib, ang_moment;
  Color col;
  float angle;
  void update(float dt)
  {
    pos += vel * dt;
    // float dt = dt_ms;
  }
  void draw(Graphics &g, Mesh &m)
  {
    g.pushMatrix();
    g.translate(pos);
    //    g.rotate(angle, vel);
    //    g.scale(vib);
    g.color(col);
    g.draw(earth_sphere);
    g.popMatrix();
  }
};

struct Spikes
{
  Vec3f pos, vib;
  Color col;
  Spikes()
  {
    pos = randomVec3f(10);
  }
  // int cov_id;
  void draw(Graphics &g, Mesh &m)
  {
    m.vertex(pos);
    m.color(col);
  }
};

struct Node
{
  Vec3f position;
  Color col;
  Node()
  {
    position = randomVec3f(70);
    //   vel =  Vec3f( 0, (rnd::uniform()-0.5) , 0);
    //col = RGB(1,1,1);
  }
  void draw(Graphics &g, Mesh &m)
  {
    m.vertex(position);
    m.color(col);
  }
};

struct Float3
{
  float data[3] = {};
  void set(float x, float y, float z)
  {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  float operator[](size_t i) const { return data[i]; }
  float &operator[](size_t i) { return data[i]; }
};


struct MyApp : public App
{
  Parameter time{"/time", "", 0.01, "", 0.001, 127};
  Parameter alignment{"/alignment", "", 0.125, "", 0.001, 0.3};
  //  Parameter gain{"gain", 0.01, 0, 1};

  Parameter frequency{"Frequency", 100, 0, 1000};
  Parameter modulation{"Modulation", 1., 0.1, 10.};
  Parameter interactor{"Pollution Index", 5, 0.1, 10};
  Parameter gain{"Gain", 1, -30, 30};

  // Boid boids[Nb];
  Mesh box;
  Mesh mesh;
  vector<Cells> cells;
  vector<earths> earths;
  vector<Spikes> spikes;
  vector<Node> node;

  ShaderProgram shader;
  float vid_radius = 10;
  float cell_radius = 3;
  float spike_radius = 0.3;
  float force;
  bool box_draw = false;
  Vec3f unitVector;
  Vec3f gravity;
  Vec3f init_vec[Nb];
  int sa, boundary;
  float sounda[422];
  float cell_angle;
  // outside is 422 inside is 1000 each
  float data[422][1000];
  Texture texture;
  float back_color_phase = 1;

  static const int senario = 8;
  Image imageData[senario];
  Texture tex[senario];

  gam::Sine<> carrier[Nb];
  gam::Sine<> modulator[Nb];
  gam::Saw<> interact_saw[Nb];
  Line line_saw;
  double a = 0;
  double b = 0;
  float halfSize = 0.1;

  bool freeze = false;
  float collision[Nb] = {0};
  float timer = 0;
  int frame[Nb];
  int init_behave[Nb];
  BufferObject buffer;
  vector<Float3> spike_positions;
  float theta[Nb][Np], beta[Nb][Np]; // spike degree
  osc::Recv server;
  int button;

  void PlatformSetupSize()
  {
    int total_width, total_height;
    al::sphere::getFullscreenDimension(&total_width, &total_height);
    std::cout << total_width << ", " << total_height << std::endl;
    dimensions(0, 0, total_width, total_height);
  }
  void onCreate()
  {    

    // // OSC comm
    server.open(7777, "0.0.0.0", 0.05);
    server.handler(oscDomain()->handler());
    server.start();
    boundary = 600;
    node.resize(nodeCount);
    spikes.resize(Nb * Np);
    // CSVReader reader;
    // reader.addType(CSVReader::REAL);
    // reader.readFile("data/Y_voltage_force_flatten_transpose.csv");
    // std::vector<double> column0 = reader.getColumn(0);
    for (int i = 0; i < 422; i++)
    {
      for (int j = 0; j < 1000; j++)
      {
        data[i][j] = 0;
      }
    }

    cell_angle = al::rnd::uniform();
    nav().pos(0, 0, -40);
    nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphere(cell_sphere, cell_radius);
    addSphereWithTexcoords(earth_sphere, vid_radius,200);
    back_mesh.primitive(Mesh::POINTS);
    spike_mesh.primitive(Mesh::POINTS);
    shader.compile(vertex1, fragment1, geometry1);

    imageData[0] = Image("../texture/EarthTexture.png");
    imageData[1] = Image("../texture/EarthTexture.png");
    imageData[2] = Image("../texture/EarthTexture.png");
    imageData[3] = Image("../texture/EarthTexture.png");
    imageData[4] = Image("../texture/EarthTexture.png");
    imageData[5] = Image("../texture/EarthTexture.png");
    imageData[6] = Image("../texture/EarthTexture.png");
    imageData[7] = Image("../texture/EarthTexture.png");
    for (int i = 0 ; i < senario; i++){
      tex[i].create2D(imageData[i].width(), imageData[i].height());
      tex[i].submit(imageData[i].array().data(), GL_RGBA, GL_UNSIGNED_BYTE);  
    }


    
    // shader.compile(slurp("../shd/point-vertex.glsl"),
    //                     slurp("../shd/point-fragment.glsl"),
    //                     slurp("../shd/point-geometry.glsl"));
    // Create Spike vertex
    for (int i = 0; i < Nb; i++)
    {
      for (int j = 0; j < Np; j++)
      {
        theta[i][j] = al::rnd::uniform(2*3.141592);
        beta[i][j] = al::rnd::uniform(2*3.141592);
      }
    }
    for (int i = 0; i < Nb; ++i)
    {
      init_vec[i] = randomVec3f(1.0);
    }

    for (int i = 0; i < nodeCount; i++)
    {
      Node &n = node[i];
      n.col = HSV(0.9, 0.0, 1);
      //      n.position = Vec3f((al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale);
    }
    for (int i = 0; i < Nb * Np; i++)
    {
      Spikes &s = spikes[i];
      s.col = HSV(1, 1, 1.0);
      //      n.position = Vec3f((al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale);
    }

    cell_sphere.generateNormals();
    earth_sphere.generateNormals();

    cells.resize(Nb);
    earths.resize(Nb);

    texture.create2D(300, 300, Texture::R8, Texture::RED, Texture::SHORT);
    int Nx = texture.width();
    int Ny = texture.height();
    std::vector<short> alpha;
    alpha.resize(Nx * Ny);
    for (int j = 0; j < Ny; ++j)
    {
      float y = float(j) / (Ny - 1) * 2 - 1;
      for (int i = 0; i < Nx; ++i)
      {
        float x = float(i) / (Nx - 1) * 2 - 1;
        float m = exp(-13 * (x * x + y * y));
        m *= pow(2, 15) - 1; // scale by the largest positive short int
        alpha[j * Nx + i] = m;
      }
    }
    texture.submit(&alpha[0]);

    resetCells(); // reset
    // Randomly distribute the starting frames
    for (int i = 0; i < Nb; i++)
    {
      // frame[i] = (int)al::rnd::uniform(999);
      frame[i] = (int)al::rnd::uniform(999);
    }
  }

  // Randomize boid positions/velocities uniformly inside unit disc
  void resetCells()
  {
    for (auto &c : cells)
    {
      c.pos = 0;
      c.vel = 0;
    }
  }

  void onAnimate(double dt)
  {

    // move nav
    Vec3f point_you_want_to_see = Vec3f(0,0,0); // examplary point that you want to see
    nav().faceToward(point_you_want_to_see, Vec3f(0, 1, 0), 0.7);
    nav().moveF( (OSCvibDepth - 64) * 0.001);
    nav().moveR( 0.1 * cos( (tableH - 64)*0.1 ));

    a += 0.001;
    b += 0.001;
    if (freeze)
      return;

    dt = time;
    timer += dt;
    //     cout << dt << endl;
    // Compute boid-boid interactions
    // vector<Vec3f> &position(mesh.vertices());
    for (int i = 0; i < Nb; ++i)
    {
      for (int j = i + 1; j < Nb; ++j)
      {
        // printf("checking cells %d and %d\n", i,j);

        auto ds = cells[j].pos - cells[i].pos;
        auto dist = ds.mag();

        // dist = position[j] - position[i];

        // Collision avoidance
        // float pushRadius = 0.05;
        // float pushStrength = 1;
        float push = 0;

        auto pushVector = ds.normalized() * push;
        cells[i].vel += pushVector;
        cells[j].vel -= pushVector;

        // Velocity matching
        // float matchRadius = 0.125;
        float nearness = exp(-al::pow2(dist / alignment));
        Vec3f veli = cells[i].vel;
        Vec3f velj = cells[j].vel;

        // Take a weighted average of velocities according to nearness
        cells[i].vel = veli * (1 - 0.5 * nearness) + velj * (0.5 * nearness);
        cells[j].vel = velj * (1 - 0.5 * nearness) + veli * (0.5 * nearness);
        // cells[i].ang_moment = cells[j].vel.mag()*10;
        cells[i].angle += cells[j].vel.mag() * 100;
        cells[i].vib = Vec3f{0.7, 0.7, 0.7} + cells[j].vel * 0.01;
      }
    }

    for (int j = 0; j < Nb; j++)
    {
      // find the center of the flock
      const Vec3f &point = cells[j].pos;

      int count = 0;
      Vec3f sum(0, 0, 0);
      cells[j].acc = {0, 0, 0};
      // for (int i = 0; i < Nb; i++)
      // {
      //   if ((point - cells[i].pos).mag() < separation)
      //   {
      //     sum += cells[i].pos; // find the center
      //     count++;
      //     // find the bounds of these points
      //     // cout << count << endl;
      //   }
      // }
      // if (count > 0)
      // {
      //   Vec3f center = sum / count; // what if count == 0?
      //   auto dsc = cells[j].pos - center;
      //   if (dsc.magSqr() > 0)
      //   {
      //     force = centering / dsc.magSqr();
      //     if (force > 1)
      //     {
      //       force = 1;
      //     }
      //     unitVector = dsc.normalized();
      //     gravity = unitVector * force;
      //     cells[j].acc += gravity;
      //     cells[j].vel += cells[j].acc;
      //   }
      // }
    }


    for (int i = 0; i < Nb; i++)
    {
      cells[i].update(dt);
      // heads.vertex(cells[i].pos);
      // heads.color(HSV(float(i) / Nb * 0.3 + 0.3, 0.7));
      cells[i].col = HSV(float(i) / Nb * 0.01 + 1, 0.5, 0.3);
      // heads.color(HSV(rnd::uniform(), 1.0f, 1.0f));

      // Define earths position
      earths[i].pos = cells[i].pos;
      earths[i].col = HSV(float(i) / Nb * 0.01 + 0.1, 0.1, 0.8);
      // detect collision
      // if( (earths[i].pos - cells[i].pos).mag() < (cell_radius - vid_radius+0.9) ){
      if( data[init_behave[i]][frame[i]] < 3 ){
        collision[i] = 1;
      }else{
        collision[i] = 0;
      }

    }
    float spike_dist = 1.01 * vid_radius;
    // Spike positions
    
    if (button != interactor){
      for (int i = 0; i < Nb; i++)
      {
        for (int j = 0; j < Np; j++)
        {
          theta[i][j] = al::rnd::uniform(2*3.141592);
          beta[i][j] = al::rnd::uniform(2*3.141592);

          spikes[Nb * (j) + i].pos = earths[i].pos + Vec3f(spike_dist * cos(theta[i][j]) * sin(beta[i][j]), spike_dist * sin(theta[i][j]) * sin(beta[i][j]), spike_dist * cos(beta[i][j]));
          spikes[Nb * (j) + i].col = HSV(1, 1, 1);
        }
      }
      button = interactor;  
    }
    time.set(OSCcarMul);
    alignment.set(OSCmodMul);
    frequency.set(OSCvibRate);
    modulation.set(OSCvibDepth);
    interactor.set(int(OSCtable)+1);
    // gui.add(time);
    // gui.add(alignment);
    // gui.add(frequency);
    // gui.add(modulation);
    // gui.add(interactor);
    // gui.add(gain);

    // // Loop the 0~1000 Frames
    // for (int i = 0; i < Nb; i++)
    // {
    //   frame[i] += 1;
    //   if (frame[i] == 999)
    //   {
    //     frame[i] = 0;
    //     init_behave[i] = (int)al::rnd::uniform(421);
    //   }
    // }
  }

  void onDraw(Graphics &g)
  {
    g.clear(0);

    g.lighting(true);
    g.depthTesting(true);
    g.blending(true);
    g.blendTrans();
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
    g.meshColor();

    // Cells
    //    cell_sphere.reset();
    {
      g.pushMatrix();
      for (auto c : cells)
      {
        // c.draw(g, cell_sphere);
      }
      g.popMatrix();
    }

    // earths
    //    earth_sphere.reset();
      g.texture();
      g.pushMatrix();
      for (auto v : earths)
      {         
        tex[0].bind();
        v.draw(g, earth_sphere);
        tex[0].unbind();
      }
      g.popMatrix();

    spike_mesh.reset(); /////////////////////////////////////////////////////////////
    // for (int i = 0; i < Nb * Np; i++)
    {
      texture.bind();
      for (auto s : spikes)
        s.draw(g, spike_mesh);
    }
    g.meshColor();
    g.shader(shader);
    g.shader().uniform("halfSize", spike_radius);
    g.draw(spike_mesh);
    texture.unbind();
    // back_mesh.reset(); /////////////////////////////////////////////////////////////
    // for (int i = 0; i < nodeCount; i++)
    {
      texture.bind();
      for (auto n : node)
        n.draw(g, back_mesh);
    }
    g.meshColor();
    g.shader(shader);
    g.shader().uniform("halfSize", halfSize);
    g.draw(back_mesh);
    texture.unbind();
  }
  float sound_fin = 0;
  float saw;

  void onSound(AudioIOData &io) override
  {
    while (io())
    {
      for (int vid = 0; vid < Nb; vid++)
      {
        carrier[vid].freq(frequency.get() + (data[init_behave[vid]][frame[vid]] * modulation.get()));
        float v = carrier[vid]() * dbtoa(gain.get());

        interact_saw[vid].freq(200 + ((50 + data[init_behave[vid]][frame[vid]]) * interactor.get()));
        line_saw.set(interact_saw[vid]() * collision[vid]);
        saw = interact_saw[vid]();
        sounda[vid] = ( v + saw )* 0.0001;
      }
      for (int vid = 0; vid < Nb; vid++)
      {
        sound_fin += sounda[vid];
      }
      io.out(0) = io.out(1) = sound_fin;
    }
  }

  bool onKeyDown(const Keyboard &k)
  {
    switch (k.key())
    {
    case 'r':
      resetCells();
      break;
    case ']':
      box_draw = !box_draw;
      break;
    case '.':
      boundary++;
      cout << "boundary: " << boundary << endl;
      break;
    case ',':
      boundary--;
      cout << "boundary: " << boundary << endl;
      break;
    case ' ':
      freeze = !freeze;
      break;
    }

    return true;
  }
  void onMessage(osc::Message &m) override
  {
    // m.print();
    if (m.addressPattern() == string("/car"))
    {
      m >> OSCcarMul;
    }
    else if (m.addressPattern() == string("/mod"))
    {
      m >> OSCmodMul;
    }
    else if (m.addressPattern() == string("/vib"))
    {
      m >> OSCvibRate;
    }
    else if (m.addressPattern() == string("/vibDepth"))
    {
      m >> OSCvibDepth;
    }
    else if (m.addressPattern() == string("/tableH"))
    {
      m >> tableH;
    }
    else if (m.addressPattern() == string("/table"))
    {
      m >> OSCtable;
    }
    else{
      m.print();
    }
  }

  void onInit() override
  {
    if (al_get_hostname() == "moxi" || fullscreen)
    {
      PlatformSetupSize();
    }
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(time);
    gui.add(alignment);
    gui.add(frequency);
    gui.add(modulation);
    gui.add(interactor);
    gui.add(gain);
   }
};

int main() { MyApp().start(); }

string slurp(string fileName)
{
  fstream file(fileName);
  string returnValue = "";
  while (file.good())
  {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}
Vec3f randomVec3f(float scale)
{
  return Vec3f(al::rnd::uniformS(), al::rnd::uniformS(), al::rnd::uniformS()) * scale;
}

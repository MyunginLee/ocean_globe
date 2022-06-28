// earth . 
#include "headers_10.hpp"

using namespace al;
using namespace std;
using namespace gam;

Mesh earth_sphere, back_mesh, point_mesh;
Pan<> mPan;
Env<3> mAmpEnv;
float OSCcarMul, OSCmodMul, OSCvibRate, OSCvibDepth, OSCtable,tableH; 

LFO<> mod;
OnePole<> onePole;
string slurp(string fileName);
Vec3f randomVec3f(float scale);

static bool fullscreen = false;

struct Points
{
  Vec3f pos, vib;
  Color col;
  Points()
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
  Parameter interactor{"Pollution Index", 0.1, 0, 8};
  Parameter gain{"Light", 1, 0, 1};

  // Boid boids[Nb];
  Mesh box;
  Mesh mesh;
  vector<Points> points;
  vector<Node> node;

  ShaderProgram shader;
  float vid_radius = 5;
  float point_radius = 0.1;
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
  Light light;          
  float light_phase = 0;
  bool freeze = false;
  float collision[Nb] = {0};
  float timer = 0;
  int frame[Nb];
  int init_behave[Nb];
  BufferObject buffer;
  vector<Float3> point_positions;
  float theta[Nb][Np], beta[Nb][Np]; // point degree
  osc::Recv server;
  int button;
  string fileName;

  Scene *ascene{nullptr};
  Vec3f scene_min, scene_max, scene_center;
  vector<Mesh> mMesh;
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
    points.resize(Nb * Np);
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
    nav().pos(0, 0, 40);
    // nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphereWithTexcoords(earth_sphere, vid_radius, 1024, 0);
    back_mesh.primitive(Mesh::POINTS);
    point_mesh.primitive(Mesh::POINTS);
    shader.compile(vertex1, fragment1, geometry1);




    imageData[0] = Image("../texture/EarthTexture.jpeg");
    imageData[1] = Image("../texture/EarthTexture2.jpeg");
    imageData[2] = Image("../texture/EarthTexture3.png");
    imageData[3] = Image("../texture/EarthTexture4.jpg");
    imageData[4] = Image("../texture/EarthTexture.jpeg");
    imageData[5] = Image("../texture/EarthTexture2.jpeg");
    imageData[6] = Image("../texture/EarthTexture3.png");
    imageData[7] = Image("../texture/direct_human_2013_impact.tif");
 

    for (int i = 0 ; i < senario; i++){
      tex[i].create2D(imageData[i].width(), imageData[i].height());
      tex[i].submit(imageData[i].array().data(), GL_RGBA, GL_UNSIGNED_BYTE);  
    }

    
    // shader.compile(slurp("../shd/point-vertex.glsl"),
    //                     slurp("../shd/point-fragment.glsl"),
    //                     slurp("../shd/point-geometry.glsl"));
    // Create point vertex
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
      Points &s = points[i];
      s.col = HSV(1, 1, 1.0);
      //      n.position = Vec3f((al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale, (al::rnd::uniform() - 0.5) * space_scale);
    }
    earth_sphere.generateNormals();

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

    // Randomly distribute the starting frames
    for (int i = 0; i < Nb; i++)
    {
      // frame[i] = (int)al::rnd::uniform(999);
      frame[i] = (int)al::rnd::uniform(999);
    }
  }

  void onAnimate(double dt)
  {

    // move nav
    Vec3f point_you_want_to_see = Vec3f(0,0,0); // examplary point that you want to see
    nav().faceToward(point_you_want_to_see, Vec3f(0, 1, 0), 0.7);
    // nav().moveF( (OSCvibDepth - 64) * 0.001);
    // nav().moveR( (tableH - 64)*0.001 );
    // nav().spinR( (tableH - 64)*0.01 );

    a += 0.001;
    b += 0.001;
    if (freeze)
      return;

    dt = time;
    timer += dt;
   
    float point_dist = 1.01 * vid_radius;
    // point positions
    
    if (button != interactor){
      for (int i = 0; i < Nb; i++)
      {
        for (int j = 0; j < Np; j++)
        {
          theta[i][j] = al::rnd::uniform(2*3.141592);
          beta[i][j] = al::rnd::uniform(2*3.141592);

          points[Nb * (j) + i].pos = Vec3f(point_dist * cos(theta[i][j]) * sin(beta[i][j]), point_dist * sin(theta[i][j]) * sin(beta[i][j]), point_dist * cos(beta[i][j]));
          points[Nb * (j) + i].col = HSV(1, 1, 1);
        }
      }
      button = interactor;  
    }
    interactor.set(int(interactor));
    OSCtable = interactor;
    // time.set(OSCcarMul);
    // alignment.set(OSCmodMul);
    // frequency.set(OSCvibRate);
    // modulation.set(OSCvibDepth);
    // interactor.set(int(OSCtable)+1);
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
    // Set light position
    light.pos(nav().pos().x,nav().pos().y, nav().pos().z);
    Light::globalAmbient({gain, gain, gain});

  }

  void onDraw(Graphics &g)
  {
    g.clear(0);

    g.lighting(true);
    g.light(light);
    g.depthTesting(true);
    g.blending(true);
    g.blendTrans();
    g.rotate(a, Vec3f(0, 1, 0));
    g.rotate(b, Vec3f(1));
    g.meshColor();
    g.texture(); // use texture

    g.translate(0,0,0);
    g.pushMatrix();
    // g.scale(3.12);
    // g.rotate(180);

    tex[int(OSCtable)].bind();
    g.draw(earth_sphere);      
    tex[int(OSCtable)].unbind();

    // for (auto &m : mMesh)
    // {
    //   tex[int(OSCtable)].bind();
    //   g.draw(m);
    //   tex[int(OSCtable)].unbind();
    // }
    g.popMatrix();

    point_mesh.reset(); /////////////////////////////////////////////////////////////
    // for (int i = 0; i < Nb * Np; i++)
    {
      texture.bind();
      for (auto s : points)
        s.draw(g, point_mesh);
    }
    g.meshColor();
    g.shader(shader);
    g.shader().uniform("halfSize", point_radius);
    g.draw(point_mesh);
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
      // resets();
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

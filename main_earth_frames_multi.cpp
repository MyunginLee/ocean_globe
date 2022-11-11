// earth . 
#include "headers_10.hpp"
#include <typeinfo>

using namespace al;
using namespace std;
using namespace gam;

Mesh earth_sphere, back_mesh, extra, target;
Mesh pic[12][2];
Vec3f img_pos, extra_pos;
int key=1;
int prevkey=0; // data to data bridge

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
  Parameter time{"Time", "", 0.01, "", 0.001, 127};
  Parameter alignment{"Alignment", "", 0.125, "", 0.001, 0.3};
  //  Parameter gain{"gain", 0.01, 0, 1};

  Parameter frequency{"Density", 100, 0, 1000};
  Parameter modulation{"Stressors ", 1., 0.1, 10.};
  Parameter interactor{"Face", 0.1, 0, 12};
  Parameter gain{"Light", 1, 0, 1};
  Parameter shift{"Shift", 0, 0, 13};
  Parameter rot{"Rotate", 1, 0, 2};
  Parameter year{"Year", 2003, 2003, 2013};


  // Boid boids[Nb];
  Mesh box;
  Mesh mesh;
  vector<Points> points;
  vector<Node> node;

  ShaderProgram shader;
  float vid_radius = 5;
  float point_radius = 0.3;
  float force;
  bool box_draw = true;
  Vec3f unitVector;
  Vec3f gravity;
  Vec3f init_vec[Nb];
  int sa;
  float cell_angle;
  // outside is 422 inside is 1000 each
  Texture texture;
  float back_color_phase = 1;
  float point_dist = 1.01 * vid_radius;
  bool draw_data = false;
  bool molph = false;
  static const int senario = 3;
  static const int years = 11;

  Image oceanData[years][2];
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
  BufferObject buffer;
  vector<Float3> point_positions;
  float theta[Nb][Np], beta[Nb][Np]; // point degree
  osc::Recv server;
  int button;
  string fileName;
  int data_W, data_H;
  Scene *ascene{nullptr};
  Vec3f scene_min, scene_max, scene_center;
  vector<Mesh> mMesh;
  Texture texBlur;
  void PlatformSetupSize()
  {
    int total_width, total_height;
    al::sphere::getFullscreenDimension(&total_width, &total_height);
    std::cout << total_width << ", " << total_height << std::endl;
    dimensions(0, 0, total_width, total_height);
  }
  void onCreate()
  {    
    // Bring ocean data (image)
    // SST
    for (int d = 0; d< years; d++){
      ostringstream ostr;
      ostr << "../texture/sst/sst_05_" << d+2003 << "_equi.png";  // ** change stressor
      char *filename = new char[ostr.str().length()+1];
      strcpy(filename, ostr.str().c_str());
      oceanData[d][0] = Image(filename);
      if (oceanData[d][0].array().size() == 0) {
        cout << "failed to load image" << endl;
        exit(1);
      }
      cout << "loaded image size: " << oceanData[d][0].width() << ", "
          << oceanData[d][0].height() << endl;
      pic[d].primitive(Mesh::POINTS);
      // pic[d].update();
    }
    // Nutrients
    for (int d = 0; d< years; d++){
      ostringstream ostr;
      ostr << "../texture/nutrient/nutrient_pollution_impact_010_" << d+2003 << "_equi.png";  // ** change stressor
      char *filename = new char[ostr.str().length()+1];
      strcpy(filename, ostr.str().c_str());
      oceanData[d][1] = Image(filename);
      if (oceanData[d][1].array().size() == 0) {
        cout << "failed to load image" << endl;
        exit(1);
      }
      cout << "loaded image size: " << oceanData[d][1].width() << ", "
          << oceanData[d][1].height() << endl;
      pic[d].primitive(Mesh::POINTS);
      // pic[d].update();
    }
    data_W = oceanData[0][1].width();
    data_H = oceanData[0][1].height();
    extra.primitive(Mesh::LINES);

    for (int d = 0; d< years; d++){
        for (int row = 0; row < data_H; row++) {
          double theta = row * M_PI /data_H;
          double sinTheta = sin(theta);
          double cosTheta = cos(theta);
          for (int column = 0; column < data_W; column++) {
            auto pixel = oceanData[d].at(column, data_H - row - 1);       
            // if( pixel.r > 80)
            //  && pixel.r < 170)
            // if( pixel.r < 170)
            if( pixel.r)
            {
            // {  
              double phi = column * M_2PI / data_W;
              double sinPhi = sin(phi);
              double cosPhi = cos(phi);

              double x = sinPhi * sinTheta;
              double y = -cosTheta;
              double z = cosPhi * sinTheta;

              double u = 1.0 - ((double)column / data_W);
              double v = (double)row / data_W;
              // target.vertex(x*point_dist,y*point_dist,z*point_dist);
              pic[d].vertex(x*point_dist + shift,y*point_dist,z*point_dist);
              // pic[d].color(HSV(pixel.r/ 100, 1, 1));
              // pic[d].color(HSV( 1 - pixel.r/ 100, 0.8-pixel.r/ 30, 0.6 + atan(pixel.r/ 240)));

              // pic[d].color(HSV( 0.6 * (pixel.r/ 100), 0.8-pixel.r/ 300, 0.6 + atan(pixel.r/ 300)));
              // pic[d].color(RGB( pixel.r/ 30, pixel.r/ 40, 0.1 + pixel.r/50 ));

              pic[d].color(RGB( 0.2 + pixel.r/30, 0.2 + pixel.r/40, 0.5 + pixel.r/50));
              // pic[d].color(HSV( 0.6 * sin(pixel.r/ 100), 0.8-pixel.r/ 300, 0.6 + atan(pixel.r/ 300)));
              // pic[d].color(HSV( pixel.r/10, 0.6 + pixel.r / 240, 0.6 + pixel.r / 240));
              // cout << pixel.r << atan(pixel.r) << endl;
          }
        }
      }
    }

    // // OSC comm
    server.open(7777, "0.0.0.0", 0.05);
    server.handler(oscDomain()->handler());
    server.start();
    node.resize(nodeCount);
    points.resize(Nb * Np);
    cell_angle = al::rnd::uniform();
    nav().pos(30, 0, 40);
    // nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphereWithTexcoords(earth_sphere, vid_radius, 1024, false);
    // earth_sphere.update();
    back_mesh.primitive(Mesh::POINTS);
    // back_mesh.update();
    shader.compile(vertex1, fragment1, geometry1);

    imageData[0] = Image("../texture/EarthTexture3.png");
    imageData[1] = Image("../texture/EarthTexture2.png");
    imageData[2] = Image("../texture/EarthTexture1.jpeg");

    for (int i = 0 ; i < senario; i++){
      tex[i].create2D(imageData[i].width(), imageData[i].height());
      tex[i].submit(imageData[i].array().data(), GL_RGBA, GL_UNSIGNED_BYTE);  
    }

    for (int i = 0; i < nodeCount; i++)
    {
      Node &n = node[i];
      n.col = HSV(0.9, 0.0, 1);
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
    texBlur.filter(Texture::LINEAR);
  }

  void onAnimate(double dt)
  {

    // move nav
    Vec3f point_you_want_to_see = Vec3f(0,0,0); // examplary point that you want to see
    nav().faceToward(point_you_want_to_see, Vec3f(0, 1, 0), 0.7);
    // nav().moveF( (OSCvibDepth - 64) * 0.001);
    // nav().moveR( (tableH - 64)*0.001 );
    // nav().spinR( (tableH - 64)*0.01 );
    // nav().moveF( (10) * 0.001);
    // nav().moveR( (10)*0.001 );



    a += 0.001;
    b += 0.001;
    if (freeze)
      return;

    dt = time;
    timer += 10*dt;
    
    
    if(molph){
      year = year + 10 *dt;
      if(year > 2013){
        year = 2003;
      }
    }
//     for (int lat = 0; lat < extra.vertices().size(); lat++) {
//       target.vertices()[lat].lerp(extra.vertices()[lat], 0.1);
//     }


    // OSC control parameters. not used for normal demo
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

    if (year && box_draw)
    {
      // if (target.vertices()[0] != extra.vertices()[0])
      {
        prevkey = key;
        key = year-2002;
        extra = pic[key-1];
        if (prevkey == key){
          target = pic[key-1];  
        } else{
          target = pic[prevkey];
        }
        year.set(year); 
      }
    }

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
    // g.rotate(80, Vec3f(0, 1, 0));

    g.translate(0,0,0);
    g.pushMatrix();
    // g.scale(3.12);

    tex[int(OSCtable)].bind();
    g.draw(earth_sphere);      
    tex[int(OSCtable)].unbind();

    g.popMatrix();

    texture.bind();
    {
      for (auto n : node)
        n.draw(g, back_mesh);
    }
    g.meshColor();
    g.shader(shader);
    g.shader().uniform("halfSize", halfSize);
    g.draw(back_mesh);
    texture.unbind();

    if(draw_data){
    // draw ocean data samples
      g.pushMatrix();
      g.blending(true);
      g.blendTrans();
      g.meshColor();
      float pointssize = nav().pos().mag();
      if (pointssize < 1)
        pointssize = 1; 
      g.pointSize(700/nav().pos().magSqr());
      // cout << 200/nav().pos().magSqr() << endl;
      g.translate(Vec3f(shift,0,0));
      g.draw(target);
      g.popMatrix();
    }    

    // texBlur.resize(fbWidth(), fbHeight());
    // g.tint(0.999);
    // g.quadViewport(texBlur, -1, -1, 2, 2);              
    // g.tint(1);  // set tint back to 1
    // texBlur.copyFrameBuffer();

  }

  void onSound(AudioIOData &io) override
  {
    while (io())
    {

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
    case ' ':
      freeze = !freeze;
      break;

//
      case '9':
        molph = !molph;
        year = 2003;
        break;

      case '0':
        draw_data = !draw_data;
        break;

        default:
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
    gui.add(shift);
    gui.add(rot);
    gui.add(year);

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

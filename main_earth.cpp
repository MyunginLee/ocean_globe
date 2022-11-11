// earth . 
#include "headers_10.hpp"
#include <typeinfo>

using namespace al;
using namespace std;
using namespace gam;

Mesh earth_sphere, back_mesh, point_mesh, extra, target;
Mesh pic[11];
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
  Parameter rot{"Scenario", 1, 0, 7};
  Parameter year{"Year", 2003, 2003, 2012};


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
  int sa, boundary;
  float sounda[422];
  float cell_angle;
  // outside is 422 inside is 1000 each
  float data[422][1000];
  Texture texture;
  float back_color_phase = 1;
  float point_dist = 1.01 * vid_radius;
  bool draw_data = false;
  bool molph = false;
  static const int senario = 3;
  static const int years = 11;

  Image oceanData[years];
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
  int data_W, data_H;
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
    // Bring ocean data (image)
    for (int d = 0; d< years; d++){
      ostringstream ostr;
      // ostr << "../texture/sst/equi_sst_" << d+2003 << ".tiff";
      // ostr << "../texture/sst/test_median/sst_" << d+2003 << ".jpg";  // ** change stressor
      // ostr << "../texture/direct_human/test_median/np_" << d+2003 << ".jpg";
      // ostr << "../texture/sst/test_median/test_0.5_sst_" << d+2003 << ".png";  // ** change stressor
      // ostr << "../texture/sst/sst_05_" << d+2003 << "_equi.png";  // ** change stressor
      ostr << "../texture/nutrient/nutrient_pollution_impact_010_" << d+2003 << "_equi.png";  // ** change stressor
      char *filename = new char[ostr.str().length()+1];
      strcpy(filename, ostr.str().c_str());
      oceanData[d] = Image(filename);
      if (oceanData[d].array().size() == 0) {
        cout << "failed to load image" << endl;
        exit(1);
      }
      cout << "loaded image size: " << oceanData[d].width() << ", "
          << oceanData[d].height() << endl;
      pic[d].primitive(Mesh::POINTS);
    }
    data_W = oceanData[0].width();
    data_H = oceanData[0].height();
    extra.primitive(Mesh::LINES);

    // // grant ocean data to pixels - Equirect
    // for (int d = 0; d< years; d++){
    //     for (int row = 0; row < data_H; row++) {
    //       double theta = row * M_PI /data_H;
    //       double sinTheta = sin(theta);pic


    //       for (int column = 0; column < data_W; column++) {
    //         auto pixel = oceanData[d].at(column, data_H - row - 1);       
          
          
    //         if( pixel.r > 0 && pixel.r < 200){
    //         // {  
    //           double phi = column * M_2PI / data_W;
    //           double sinPhi = sin(phi);
    //           double cosPhi = cos(phi);

    //           double x = sinPhi * sinTheta;
    //           double y = cosTheta * cosTheta;
    //           double z = cosPhi * sinTheta;

    //           double u = 1.0 - ((double)column / data_W);
    //           double v = (double)row / data_W;
    //           // target.vertex(x*point_dist,y*point_dist,z*point_dist);
    //           pic[d].vertex(x*point_dist + shift,y*point_dist,z*point_dist);
    //           pic[d].color(HSV( 0.6 + pixel.r/ 240, 0.6 + atan(pixel.r/ 240), 0.6 + atan(pixel.r/ 240)));

    //           // pic[d].color(HSV( atan(pixel.r/30), 0.6 + pixel.r / 240, 0.6 + pixel.r / 240));
    //           // cout << pixel.r << atan(pixel.r) << endl;
    //       }
    //     }
    //   }
    // }


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
            if( pixel.r > 0)
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
              pic[d].color(RGB( pixel.r/ 10, pixel.r/ 20, pixel.r/30 ));

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
    nav().pos(30, 0, 40);
    // nav().faceToward(Vec3d(0, 0, 0), Vec3d(0, 1, 0));
    addSphereWithTexcoords(earth_sphere, vid_radius, 1024, false);
    back_mesh.primitive(Mesh::POINTS);
    point_mesh.primitive(Mesh::POINTS);
    shader.compile(vertex1, fragment1, geometry1);

    imageData[0] = Image("../texture/EarthTexture3.png");
    imageData[1] = Image("../texture/EarthTexture2.png");
    imageData[2] = Image("../texture/EarthTexture1.jpeg");

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
        theta[i][j] = al::rnd::uniform(2*M_PI);
        beta[i][j] = al::rnd::uniform(2*M_PI);
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
    // nav().moveF( (10) * 0.001);
    // nav().moveR( (10)*0.001 );



    a += 0.001;
    b += 0.001;
    if (freeze)
      return;

    dt = time;
    timer += 200*dt;
    
    
    if(molph){
      year = year + dt;
      if(year > 2013){
        year = 2003;
      }
    }
//     for (int lat = 0; lat < extra.vertices().size(); lat++) {
//       target.vertices()[lat].lerp(extra.vertices()[lat], 0.1);
//     }

// //  Adjust color
//     // point positions
    
//     if (button != interactor){
//       for (int i = 0; i < Nb; i++)
//       {
//         for (int j = 0; j < Np; j++)
//         {
//           theta[i][j] = al::rnd::uniform(2*M_PI);
//           beta[i][j] = al::rnd::uniform(2*M_PI);

//           points[Nb * (j) + i].pos = Vec3f(point_dist * cos(theta[i][j]) * sin(beta[i][j]), point_dist * sin(theta[i][j]) * sin(beta[i][j]), point_dist * cos(beta[i][j]));
//           points[Nb * (j) + i].col = HSV(1, 1, 1);
//         }
//       }
//       button = interactor;  
//     }
//     interactor.set(int(interactor));
//     OSCtable = interactor;

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

    // if (year && box_draw)
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

    point_mesh.reset(); ////////////////////////////////////////////////////////////alive 
    
    // if(box_draw){
    //   texture.bind();
    //   for (auto s : points)
    //     s.draw(g, point_mesh);
      
    //   g.meshColor();
    //   g.shader(shader);
    //   g.shader().uniform("halfSize", point_radius);
    //   g.draw(point_mesh);
    //   texture.unbind();
    // }
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


    if(draw_data){
    // draw ocean data samples
      g.pushMatrix();
      g.blending(true);
      g.blendTrans();
      g.meshColor();
      float pointssize = nav().pos().mag();
      if (pointssize < 1)
        pointssize = 1; 
      g.pointSize(500/nav().pos().magSqr());
      g.translate(Vec3f(shift,0,0));
      g.draw(target);
      g.popMatrix();
    }    
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
    // gui.add(alignment);
    gui.add(frequency);
    gui.add(modulation);
    // gui.add(interactor);
    gui.add(gain);
    // gui.add(shift);
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

#include <iostream>
#include <cmath>
#include <cstdio>
#include <allegro.h>
#include <pthread.h>
#include <fstream>
#include "Vect2.h"
using namespace std;

#define M_WIDTH 480
#define M_HEIGHT 480

int beban = 8;
double PI = 3.1415926535;
const int N = 10;
double M = 0.0152;
double L = 0.3;
double t = 0.0;
void* draw(void *arg);
void init();
void updateForce(double l, double k, double kb, double gamma, double gammab);
void updatePosition(double dt);
void writeFile();
void draw();
Vect2 cog();

pthread_t tid;

class particle
{
public:
    double radius;
    double massa;
    double teta, alfa, omega, omegadis, tetadis, Fgtan, Fgrad;
    double Tau;
    double panjang;
    Vect2 r, rt, v, a, Fs, vdis;
    particle();
    ~particle();
} par[N];

particle::particle()
{
    radius = 0.5*L/(N-1);
    massa = M/N;
    teta = 0.0;
    alfa = 0.0;
    omega = 0.0;
    omegadis = 0.0;
    tetadis = 0.0;
    Tau = 0.0;
    Fgtan = 0.0;
    Fgrad = 0.0;
    panjang = 0.0;
}

particle::~particle()
{
}

int main()
{

    double dt = 0.000001;
    double l = L/(N-1);
    double kb = 0.0058;//0.5
    double k = 4000000.0;
    double gamma = 0.001;
    double gammab = 0.0005;//0.005005
    init();

    allegro_init();
    install_keyboard();
    set_gfx_mode(GFX_AUTODETECT_WINDOWED, M_WIDTH, M_HEIGHT, 0, 0);

    //  ====== Draw THREAD ======
    int err = pthread_create(&(tid), NULL, &draw, NULL);
    if (err != 0)
        cout << "\ncan't create thread : " << (err) << endl;
    else
        cout << "\n Thread created successfully\n" << endl;
// =====================================

    int writemarker = 0;

    while(true)
    {
        if(key[KEY_ESC])
            break;

        if(key[KEY_W])
        {
            if(writemarker == 0)
                writeFile();
            writemarker = 1;
        }


        if(key[KEY_P])
        {
            while(key[KEY_P]);
        }

        updateForce(l,k,kb,gamma,gammab);
        updatePosition(dt);

        t += dt;
        //rest(1);
    }

    return 0;
}

void init()
{
    t = 0.0;
    for(int i=0; i<N; i++)
    {
        par[i].r.X = 0.5+i*(L/(N-1));
        par[i].r.Y = 0.5;
        par[i].rt = par[i].r;
    }
    par[N-1].massa *= beban;
}

void updateForce(double l, double k, double kb, double gamma, double gammab)
{
    Vect2 grav(0.0,-10.0);
    for(int i=2; i<N; i++)
    {
        Vect2 a = par[i-2].r - par[i-1].r;
        Vect2 b = par[i].r - par[i-1].r;


        par[i].Fgtan = -grav.Length() * par[i].massa * sin(grav.AngleBetween(b)*PI/180);
        par[i].Fgrad = -grav.Length() * par[i].massa * cos(grav.AngleBetween(b)*PI/180);

        par[i].teta = b.AngleBetween(a);

        double Fspringb = -kb*(par[i].teta);
        par[i].Tau = -gammab*par[i].omegadis + par[i].Fgtan*b.Length();
        par[i].Tau += Fspringb;
        par[i-1].Tau -= Fspringb;
    }
    for(int i=1; i<N-1; i++)
    {
        Vect2 c = par[i].r - par[i+1].r;
        Vect2 Fspring = (c/c.Length())*(-k*(c.Length() - l));
        par[i].Fs += Fspring;
        par[i+1].Fs = Fspring*-1.0 - par[i+1].vdis*gamma;
    }
}

void updatePosition(double dt)
{
    for(int i=2; i<N; i++)
    {
        Vect2 a = par[i].r - par[i-1].r;
        par[i].aa = a;
        par[i].a = ( (a/a.Length())*par[i].Fgrad + par[i].Fs )/par[i].massa;
        par[i].alfa = par[i].Tau/(a.Length()*a.Length()*par[i].massa);

        par[i].v = par[i].v + par[i].a*dt;
        par[i].r = par[i].r + par[i].v*dt;

        par[i].omega = par[i].omega + par[i].alfa*dt;

        int j = i;
        //for(int j=i; j<N; j++)
        {
            Vect2 b = par[j].r - par[i-1].r;
            Vect2 temp;
            temp.X = b.X*cos(par[i].omega*dt*PI/180) + b.Y*sin(par[i].omega*dt*PI/180);
            temp.Y = -b.X*sin(par[i].omega*dt*PI/180) + b.Y*cos(par[i].omega*dt*PI/180);
            par[j].r = par[i-1].r + temp;
        }
        par[i].vdis = (par[i].r-par[i].rt)/dt;

        par[i].tetadis = par[i].r.AngleBetween2(par[i].rt);
        par[i].omegadis = par[i].tetadis/dt;

        par[i].rt = par[i].r;
    }
}

void* draw(void *arg)
{
    double scale = M_WIDTH/1;
    FONT *font1 = load_font("font4.pcx", NULL, NULL);
    char time[20];
    char posY[40];

    BITMAP *buffer;
    buffer = create_bitmap(M_WIDTH,M_HEIGHT);

    while(true)
    {
        clear_to_color(buffer, makecol(255,255,255));

        for(int i = 0; i < N; i++)
        {
            double x = par[i].r.X*scale;
            double y = M_HEIGHT - par[i].r.Y*scale;
            double r = par[i].radius*scale;
            circlefill(buffer, (int)x, (int)y, (int)r, makecol(0,0,0));
        }
        sprintf(time,"t: %lf",t);
        sprintf(posY,"y: %lf",(0.5 - par[N-1].r.Y)*100);
        textout_ex(buffer, font1, time, 10, 10, makecol(0,0,0), -1);
        textout_ex(buffer, font1, posY, 10, 25, makecol(0,0,0), -1);

        acquire_screen();
        draw_sprite(screen,buffer,0,0);
        release_screen();
        //rest(1);
    }
}

void writeFile()
{
    char x[100];

    sprintf(x,"data/graph-%d.dat",beban);
    ofstream write(x);

    for(int i=0; i<N; i++)
    {
        sprintf(x,"%lf %lf\n",par[i].r.X, par[i].r.Y);
        write << x;
    }
    write.close();

    sprintf(x,"data/time.dat");
    write.open(x, ios::app | ios::out);
    sprintf(x,"%d %lf\n",beban, t);
    write << x;
    write.close();

    sprintf(x,"data/cog.dat");
    write.open(x, ios::app | ios::out);
    sprintf(x,"%d %lf %lf\n",beban, cog().X, cog().Y);
    write << x;
    write.close();
}

Vect2 cog()
{
    Vect2 sum;
    double sumMassa = 0.0;
    for(int i=0; i<N; i++)
    {
        sum += par[i].r*par[i].massa;
        sumMassa += par[i].massa;
    }
    sum /= sumMassa;

    return sum;
}

void draw()
{
    FILE *gp;
    gp=popen("gnuplot -geometry 400x400+450+0 > /dev/null 		2>&1","w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'sim-%lf.png'\n",t);
    fprintf(gp, "set xrange [0:1]\n");
    fprintf(gp, "set yrange [0:1]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set label 1 't: %lf' at 0.03,%lf\n",t,(double)1-0.03);

    for(int i=0; i<N; i++)
    {
        fprintf(gp, "set object %d circle center %lf,%lf size %lf\n", i+1, par[i].r.X, par[i].r.Y, par[i].tetadis);
    }

    fprintf(gp,"plot 0\n");
    pclose(gp);
}

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <stdlib.h>

using namespace std;
void IS (float,float, float, float, float, float &);
void GU (float, float, float, float, float &);
void Bern (float, float, float, float &);
void TG(float, float, float, float, float &);
void PV(float, float, float, float &);
void TP (float, float, float &);
void Co(float, float, float, float, double &);
void IE (float, float, float &);
void Eim(float, float, float, float &);
void Ilu (float, float, float &);
void CP(float, float, float &);
void CA(float, float, float, float, float &);

int main(int argc, char *argv[])
{
    /*empieza programa*/
    char sn[1];
    int op, opf, ope, opv;
    float f, A, K, K1, a, b, c, p, v, I, pi, m, m1, d, F, v1, p1, h, g, P, n, m2, U, V, P1, q, q1, r, X, R, Z, A2, B, mu0, F3, dF, dS, Ev, p2, r1, z, Tz, x, y, f1, E, f4, y5, Cc;
    double  F1;
    sn[0]='s';

    cout<<"De que materia quieres resolver la formula: \n";
    cout<<"\n Presiona 1 para fisica \n Presione 2 para electronica\n Presione 3 para varias\n 0 para salir "<<endl;
    cout<<"\nQue elijes?: ";
    cin>>op;
    system("cls");

    while (!((op==1)||(op==2)||(op==3)||(op==0)))
    {
        cout<<"***Ese numero no se encuentra en la lista.Intentalo de nuevo. \n"<<endl;
        cout<<"De que materia quiere resolver la formula: \n"<< endl;
        cout<<"\n Presione 1 para fisica \n presione 2 para electronica\n presione 3 para varias\n cero para salir";
        cin>>op;
        system("cls");
    }

    if(!(op==0))
    {
        switch(op)
        {
        case 1:
            cout<<"\n 1 para la formula de la intensidad sonora";
            cout<<"\n 2 para la formula de la ley de gravitacion universal";
            cout<<"\n 3 para la formula del Teorema De Bernoulli";
            cout<<"\n 4 para la formula de la Teoria Cinetica De Los Gases \n";
            cout<<"\n Que elijes?: ";
            cin>>opf;
            system("cls");
            switch(opf)
            {
            case 1:

                cout<<"La formula de la Intensidad sonora es: I=(2pi)*(f^2)*(A^2)*p*v \n";
                cout<<"\nDonde: \n f=frecuencia de vibracion \n A= amplitud de onda \n p=densidad del medio \n v=velocidad del sonido \n";

                IS (pi, v, p, A, f, I);
                cout<<"\n El valor la intensidad sonora es: " << I <<" \n";
                break;
                system("cls");
            case 2:
                cout<<"La formula de la Ley Gravitacional es: F=K(m*m1)/d^2 \n";
                cout<<"\nDonde:\nF=fuerza de atraccion\nm=masa de los cuerpos\nd=distancia entre los dos cuerpos\nK=constante  universal de la gravitacion \n";

                GU (m, m1, d, K, F);
                cout << " \n El valor de la Fuerza es: " << F <<"\n ";

                break;
                system("cls");
            case 3:
                cout<<"La formula de Bernoulli es: constante P =((v^2*p1)/2)+h*p1*g \n";
                cout<<"\nDonde:\nP=presion que se ejerce sobre el punto considerado\nv=velocidad del liquido\np1=densidad del liquido\nh=altura\ng=gravedad \n";

                Bern (v1, p1, h, P);
                cout<<" \n El valor de la Presion es: " << P << " \n";
                break;

                system("cls");
            case 4:
                cout<<"La formula para Cinetica de los gases es:P =(1/3*n*m*U^2)/V \n";
                cout<<"\nDonde:\nP=presion de gas\nn=numero de moleculas contenidas\nm=masa de la molecula\nU=velocidad de la masa\nV=volumen \n";

                TG(n, m, U, v, P1);
                cout<<"\n El valor de la Presion es:  "<< P1<< "\n";
                break;
            default:
                cout<<"opcion no valida";
            }
            break;
            system("cls");
        case 2:
            cout<<"\n 1 para la formula de la Ley De Coulomb";
            cout<<"\n 2 para la formua de la Impedancia electrica";
            cout<<"\n 3 para la formula de la Fuerza de los Electroimanes ";
            cout<<"\n 4 para la formula de la Iluminancia \n";
            cout<<"\n Que elijes?: ";
            cin>>ope;
            system("cls");
            switch(ope)
            {
            case 1:
                cout<<"La formula de la Ley de Coulomb es: F = K((q*q1)/(r^2)) \n";
                cout<<"\nDonde:\nF=fuerza entre dos cargas\nK=constante de\nq=carga 1\nq1=carga 2\nr=distancia que los separa \n";
                Co(q, q1, r, K1, F1);
                cout<<"\n El valor de la Fuerza es:  "<< F1 << "\n";
                break;
                system("cls");
            case 2:
                cout<<"La formula de Impedancia electrica es:Z =(((R^2)+(X^2))^-2 \n";
                cout<<"\nDonde: \nZ=Impedancia\nR=Resistencia\nX=Reactancia\n";
                IE(X, R, Z);
                cout<<"El valor de la Impedancia es: "<< Z << "\n";
                break;
                system("cls");
            case 3:
                cout<<"\nLa formula de la Fuerza de electroimanes es: F =(((B^2)A)/2*mu0)\n";
                cout<<"\nDonde:\nF=Fuerza\nB=Campo Magnetico\nA=Area de las Caras de los Imanes\nmu0=permeabilidad de el espacio libre \n";
                Eim(B, A2, mu0, F3);
                cout<<"\n El valor de la Fuerza es: "<< F3 << "\n";
                break;
                system("cls");
            case 4:
                cout<<"La formula de Iluminancia es: Ev = dF/dS \n";
                cout<<"\nDonde:\nEv=Iluminancia\ndF=Flujo luminoso\ndS=diferencial del area considerada \n";
                Ilu(dF, dS, Ev);
                cout<<"\n El valor de la Iluminancia es: "<< Ev << "\n";
                break;
            default:
                cout<<"opcion no valida";
            }
            break;

            system("cls");
        case 3:
            cout<<"\n 1 para la formula de la Presion Vertical";
            cout<<"\n 2 para la formula del Teorema de Pitagoras";
            cout<<"\n 3 para la formula del Cable Parabolico ";
            cout<<"\n 4 para la formula de la Carga Axial";
            cout<<"\n Que elijes?: \n";
            cin>>opv;
            system("cls");
            switch(opv)
            {
            case 1:
                cout<<"La formula de Presion Vertical es: Tz =(p/2pi)*((3(z^3))/(((r^2)+(z^2))^(5/2)) \n";
                cout<<"\nDonde\n Tz = incremento de presion vertical\nz=profundidad\nr=distancia radial \n p= presion inicial \n";

                PV (p, z, r, Tz);
                cout<<"\n El valor del incremento de Presion Vertical es: "<< Tz << "\n";
                break;
                system("cls");
            case 2:
                cout<<"\nLa formula de Teorema de Pitagoras es: c =(a^2+b^2)^2 \n";
                cout<<"\nDonde: \nc=hipotenusa\na=cateto opueto\nb=cateto adyacente \n";

                TP (a, b, c);
                cout<<"el valor de la hipotenusa es: "<< c << "\n";
                break;
                system("cls");
            case 3:
                cout<<"\nLa formula de Cable Parabolico es: f = x^2/4*y \n";
                cout<<"\nDonde:\nf=distancia focal\nx=distancia horizontal\ny=distancia vertical\n";
                CP(x, y, f1);
                cout<<"\nEl valor de la distancia focal es:  "<< f1 << "\n";
                break;
                system("cls");
            case 4:
                cout<<"\nLa formula de Carga Axial es: Cc=((2*pi^2)E)/(f*y)\n ";
                cout<<"\nDonde:\nCc=carga axial\nE=modulo de elasticidad\nf=deformacion \ny=logitud \n";
                CA(E, f4, y5, pi, Cc);
                cout<<"\nEl valor de la Carga Axial es: "<< Cc << "\n";
                break;
            default:
                cout<<"opcion invalida";

            }
            break;
        }

    }

    cout<< "\n Hasta luego!!\n";
    system ("pause");
    return 0;

}

void IS (float pi, float a, float b, float x, float y, float &c)
{
    cout<<"\ndame el valor de f: ";
    cin>>b;
    cout<<"\ndame el valor de A: ";
    cin>>x;
    cout<<"\ndame el valor de p: ";
    cin>>y;
    cout<<"\ndame el valor de v: ";
    cin>>a;
    pi=3.1416;
    c=((2*pi)*pow(b,2)*pow(x,2))*y*a;
}
void GU (float m, float m1, float d, float K, float &F)
{
    cout<<"\ndame el valor de la masa m: ";
    cin>>m;
    cout<<"dame el valor de la masa m1: ";
    cin>>m1;
    cout<<"\ndame el valor de la distancia: ";
    cin>>d;
    cout<<"\ndame el valor de K: ";
    cin>>K;
    F=K*(m*m1)/(d*d);
}

void Bern (float v1, float p1, float h, float &P)
{
    float g;

    cout<<"\ndame el valor de v: ";
    cin>>v1;
    cout<<"dame el valor de p1: ";
    cin>>p1;
    cout<<"dame el valor de h: ";
    cin>>h;
    g=9.8;
    P=((pow(v1,2)*p1)+(2*h*p1*g))/2;
}

void TG (float n, float m, float U, float v, float &P1)
{
    cout<<"\ndame el valor de n: ";
    cin>>n;
    cout<<"dame el valor de m: ";
    cin>>m;
    cout<<"dame el valor de U: ";
    cin>>U;
    cout<<"dame el valor de v: ";
    cin>>v;
   P1=(n*m*pow(U,2))/(3*v);

}

void PV (float p, float z, float r, float &Tz)
{
    float pi;
    cout<<"\ndame el valor de p: ";
    cin>>p;
    cout<<"dame el valor de z: ";
    cin>>z;
    cout<<"dame el valor de r: ";
    cin>>r;
    pi=3.1416;
    Tz=(p/2*pi)*(3*pow(z,3))/sqrt(pow((pow(r,2)+pow(z,2)),5));
}

void TP (float a, float b, float &c)
{
    cout<<"\nDame el valor de a: ";
    cin>>a;
    cout<<"\nDame el valor de b: ";
    cin>>b;
    c=sqrt((a*a)+(b*b));
}
void Co (float q, float q1, float r, float K1, double &F1)
{
    cout<<"\ndame el valor de q: ";
    cin>>q;
    cout<<"Dame el valor de q1: ";
    cin>>q1;
    cout<<"Dame el valor de r: ";
    cin>>r;
    cout<<"dame el valor de K: ";
    cin>>K1;
    F1=K1*((q*q1)/pow(r,2));
}
void IE (float X, float R, float &Z)
{
   cout<<"\ndame el valor de X: ";
                cin>>X;
                cout<<"dame el valor de R: ";
                cin>>R;
                Z=sqrt(pow(R,2)+pow(X,2));
}
void  Eim(float B, float A2, float mu0, float &F3)
{
cout<<"\n Dame el valor de B: ";
                cin>>B;
                cout<<"\nDame el valor de A: ";
                cin>>A2;
                cout<<"\nDame el valor de mu0: ";
                cin>>mu0;
                F3=((B*B)*A2)/(2*mu0);
}
void Ilu (float dF, float dS, float &Ev)
{
    cout<<"\nDame el valor de dF: ";
                cin>>dF;
                cout<<"\nDame el valor de dS: ";
                cin>>dS;
                Ev=dF/dS;
}

void CP(float x, float y, float &f1)
{
    cout<<"\ndame el valor de x: ";
                cin>>x;
                cout<<"\nDame el valor de y: ";
                cin>>y;
                f1=(x*x)/(4*y);
}
void CA(float E, float f4, float y5, float pi, float &Cc)
{
  cout<<"\nDame el valor de E: ";
                cin>>E;
                cout<<"\nDame el valor de f: ";
                cin>>f4;
                cout<<"\nDame el valor de y: ";
                cin>>y5;
                pi=3.1416;
                Cc=(((2*pi)*(2*pi)*E))/(f4*y5);
}

#include <vector>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
#include <fstream>

using namespace std;

int Nx = 100, Ny = 100;
double Lx = 1., Ly = 1., D = 1.;
int np;

double dx = Lx / (Nx + 1.);
double dy = Ly / (Ny + 1.);

//----------------------Fonction charge-----------------------------------------------

void charge(int me, int N_y, int np, int *ligne_beg, int *ligne_end)
{
    int r = N_y % np;
    int m = N_y / np;

    if (me < r)
    {
        *ligne_beg = me * (m + 1);
        *ligne_end = (me + 1) * (m + 1) - 1;
    }
    else
    {
        *ligne_beg = r + me * m;
        *ligne_end = *ligne_beg + m - 1;
    }
}

//----------Fonction qui effectue l'opération de soustraction de deux vecteurs------
vector<double> Sousvect(std::vector<double> x, std::vector<double> y)
{
    int n = y.size();
    vector<double> s(n);
    for (int i = 0; i < n; i++)
    {
        s[i] = x[i] - y[i];
    }
    return s;
}

//----------Fonction qui effectue l'opération d'addition de deux vecteurs------
vector<double> Sommevect(vector<double> x, vector<double> y)
{
    int n = x.size();
    vector<double> s(n);
    for (int i = 0; i < n; i++)
    {
        s[i] = x[i] + y[i];
    }
    return s;
}

//----------Fonction qui effectue le produit scalaire de deux vecteurs------
double produit_scalaire(vector<double> x, vector<double> y)
{
    int n = x.size();
    double ps;
    ps = 0;
    for (int i = 0; i < n; i++)
    {
        ps = ps + x[i] * y[i];
    }
    return ps;
}

//----------Fonction qui effectue le produit d'un scalaire avec un vecteur------
vector<double> multi_scalaire(double a, vector<double> x)
{
    int n = x.size();
    vector<double> s(n);
    for (int i = 0; i < n; i++)
    {
        s[i] = a * x[i];
    }
    return s;
}

//---------------------Produit matrice/vecteur----------------------------------

vector<double> Produit_matriciel_robin(int taille_recouvrement, double alpha, double beta, double gamma, double coeff_Robin, vector<double> U, int me)
{

    int ibeg, iend;
    charge(me, Ny+1, np, &ibeg, &iend);

    // domaine bas (grille)
    if (me == 0)
    {
        int N0 = (iend - ibeg + 1) + taille_recouvrement;
        int taille_0 = N0 * Nx; // la taille du vecteur en 2D du 1er domaine

        vector<double> prod(taille_0);

        // calcul du produit mat_vec
        for (int i = 0; i < taille_0; i++)
        {
            if (i >= 0 and i < Nx)
            {
                if (i == 0)
                {
                    prod[i] = alpha * U[i] + beta * U[i + 1] + gamma * U[i + Nx];
                }

                else if (i == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + alpha * U[i] + gamma * U[i + Nx];
                }

                else
                {
                    prod[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * U[i + Nx];
                }
            }

            else if (i >= Nx and i < taille_0 - Nx)
            {
                if (i % Nx == 0)
                {
                    prod[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else if (i % Nx == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else
                {
                    prod[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // else  // dernier bloque du domaine
            if (i >= taille_0 - Nx and i < taille_0)
            {
                if (i == taille_0 - Nx)
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * U[i + 1] + 2 * gamma * (U[i - Nx]);
                }

                else if (i == taille_0 - 1)
                {
                    prod[i] = beta * U[i - 1] + (alpha + coeff_Robin) * U[i] + 2 * gamma * U[i - Nx];
                }

                else
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * (U[i + 1] + U[i - 1]) + 2 * gamma * U[i - Nx];
                }
            }
        }

        return prod;
    }

    // domaine haut (grille)
    else if (me == np - 1)
    {
        int Nnp = iend - ibeg + 1;
        int taille_np = Nnp * Nx; // la taille du vecteur en 2D du domaine du haut

        vector<double> prod(taille_np);

        // calcul du produit mat_vec
        for (int i = 0; i < taille_np; i++)
        {
            if (i >= 0 and i < Nx)
            {
                if (i == 0)
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * U[i + 1] + 2 * gamma * (U[i + Nx]);
                }

                else if (i == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + (alpha + coeff_Robin) * U[i] + 2 * gamma * U[i + Nx];
                }

                else
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * (U[i + 1] + U[i - 1]) + 2 * gamma * U[i + Nx];
                }
            }

            else if (i >= Nx and i < taille_np - Nx)
            {
                if (i % Nx == 0)
                {
                    prod[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else if (i % Nx == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else
                {
                    prod[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // else  // dernier bloque du domaine
            if (i >= taille_np - Nx and i < taille_np)
            {
                if (i == taille_np - Nx)
                {
                    prod[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i - Nx]);
                }

                else if (i == taille_np - 1)
                {
                    prod[i] = beta * U[i - 1] + alpha * U[i] + gamma * U[i - Nx];
                }

                else
                {
                    prod[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * U[i - Nx];
                }
            }
        }

        return prod;
    }

    // domaine milieu
    else
    {
        int Nme = (iend - ibeg + 1) + taille_recouvrement + 1;
        int taille_me = Nme * Nx; // la taille du vecteur en 2D du kieme domaine

        vector<double> prod(taille_me);

        // calcul du produit mat_vec
        for (int i = 0; i < taille_me; i++)
        {
            if (i >= 0 and i < Nx)
            {
                if (i == 0)
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * U[i + 1] + 2 * gamma * U[i + Nx];
                }

                else if (i == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + (alpha + coeff_Robin) * U[i] + 2 * gamma * U[i + Nx];
                }

                else
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * (U[i + 1] + U[i - 1]) + 2 * gamma * U[i + Nx];
                }
            }

            else if (i >= Nx and i < taille_me - Nx)
            {
                if (i % Nx == 0)
                {
                    prod[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else if (i % Nx == Nx - 1)
                {
                    prod[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }

                else
                {
                    prod[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // else  // dernier bloque du domaine
            if (i >= taille_me - Nx and i < taille_me)
            {
                if (i == taille_me - Nx)
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * U[i + 1] + 2 * gamma * (U[i - Nx]);
                }

                else if (i == taille_me - 1)
                {
                    prod[i] = beta * U[i - 1] + (alpha + coeff_Robin) * U[i] + 2 * gamma * U[i - Nx];
                }

                else
                {
                    prod[i] = (alpha + coeff_Robin) * U[i] + beta * (U[i + 1] + U[i - 1]) + 2 * gamma * U[i - Nx];
                }
            }
        }

        return prod;
    }
}

//---------------------Fonctions f, g, h : termes source et conditions aux bords selon les cas--------------------------------

double f(double x, double y, double t, int cas)
{

    double r_f;
    if (cas == 1)
    {
        r_f = 2 * (y - pow(y, 2) + x - pow(x, 2));
    }

    else if (cas == 2)
    {
        r_f = sin(x) + cos(y);
    }

    else if (cas == 3)
    {
        r_f = exp(-pow(x - Lx / 2., 2)) * exp(-pow(y - Ly / 2., 2)) * cos(M_PI * t / 2.);
    }

    return r_f;
}

double g(double x, double y, int cas)
{
    double r_g;
    if (cas == 1)
    {
        r_g = 0.;
    }

    else if (cas == 2)
    {
        r_g = sin(x) + cos(y);
    }

    else if (cas == 3)
    {
        r_g = 0;
    }

    return r_g;
}

double h(double x, double y, int cas)
{
    double r_h;
    if (cas == 1)
    {
        r_h = 0;
    }

    else if (cas == 2)
    {
        r_h = sin(x) + cos(y);
    }

    else if (cas == 3)
    {
        r_h = 1;
    }

    return r_h;
}

//---------------------vecteur second membre----------------------------------

vector<double> second_membre(double beta, double gamma, double coeff_Robin, vector<double> xcoord, vector<double> ycoord, double t, vector<double> stencil_up, vector<double> stencil_down, int me, int taille_recouvrement, int cas)
{

    int sizex = xcoord.size();
    int sizey = ycoord.size();

    vector<double> S(sizex * sizey);

    for (int j = 1; j <= sizey; j++)
    {
        for (int i = 1; i <= sizex; i++)
        {
            S[(j - 1) * Nx + (i - 1)] = f(xcoord[i - 1], ycoord[j - 1], t, cas);
        }
    }

    int taille = sizey * sizex; // la taille du vecteur du domaine

    // cas domaine bas
    if (me == 0)
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;

            if (j == 1)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * g(xcoord[i - 1], 0., cas);
                    
                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else if (j >= 2 && j < sizey)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * (stencil_up[2 * Nx + i - 1] - stencil_up[i - 1]) + coeff_Robin * stencil_up[Nx + i - 1];

                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
        }

        return S;
    }

    // cas domaine haut
    else if (me == np - 1)
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;

            if (j == 1)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * (stencil_down[i - 1] - stencil_down[2 * Nx + i - 1]) + coeff_Robin * stencil_down[Nx + i - 1];

                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else if (j >= 2 && j < sizey)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * g(xcoord[i - 1], Ly, cas);

                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
        }

        return S;
    }

    else
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;

            if (j == 1)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * (stencil_down[i - 1] - stencil_down[2 * Nx + i - 1]) + coeff_Robin * stencil_down[Nx + i - 1];

                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else if (j >= 2 && j < sizey)
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
            else
            {
                if (i >= 1 and i < Nx + 1) // for (int i = 1; i < Nx+1 ; i++)
                {
                    S[(j - 1) * Nx + (i - 1)] += -gamma * (stencil_up[2 * Nx + i - 1] - stencil_up[i - 1]) + coeff_Robin * stencil_up[Nx + i - 1];

                    if (i == 1)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(0, ycoord[j - 1], cas);
                    }

                    else if (i == Nx)
                    {
                        S[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, ycoord[j - 1], cas);
                    }
                }
            }
        }

        return S;
    }
}


vector<double> BiCGStab(int taille_recouvrement, double beta, double alpha, double gamma, double coeff_Robin, vector <double> b, int me)
{

    int size = b.size();
    vector<double> r(size), rs(size), x(size), v(size), p(size), s(size), t(size);
    int iteration;
    int iterationmax = 10000;
    double epsilon = 1e-10;
    double rho, rho_prev;
    double alpha_GC, beta_GC, omega_GC;

    double norme_r, norme_b;

    iteration = 0;

    // initialisation des vecteurs x, v et p
    for (int i = 0; i < size; i++)
    {
        x[i] = 0.;
        v[i] = 0.;
        p[i] = 0.;
    }

    alpha_GC = 1.;
    beta_GC = 1.;
    omega_GC = 1.;

    r = Sousvect(b, Produit_matriciel_robin(taille_recouvrement, alpha, beta, gamma, coeff_Robin, x, me));
    rs = r;

    rho = 1.;

    norme_r = sqrt(produit_scalaire(r, r));
    norme_b = sqrt(produit_scalaire(b, b));

    
       while (norme_r > epsilon * norme_b && iteration < iterationmax)
        {
            rho_prev = rho;
            rho = produit_scalaire(rs, r);
            beta_GC = (rho / rho_prev) * (alpha_GC / omega_GC);

            vector<double> vec1(size), vec2(size), vec3(size);
            vec1 = multi_scalaire(omega_GC, v);
            vec2 = Sousvect(p, vec1);
            vec3 = multi_scalaire(beta_GC, vec2);
            p = Sommevect(r, vec3);

            v = Produit_matriciel_robin(taille_recouvrement, alpha, beta, gamma, coeff_Robin, p, me); // p,me
            alpha_GC = rho / produit_scalaire(rs, v);

            vector<double> vec4(size);
            vec4 = multi_scalaire(alpha_GC, v);
            s = Sousvect(r, vec4);

            t = Produit_matriciel_robin(taille_recouvrement, alpha, beta, gamma, coeff_Robin, s, me); // s,me
            omega_GC = produit_scalaire(t, s) / produit_scalaire(t, t);

            vector<double> vec5(size), vec6(size), vec7(size);
            vec5 = multi_scalaire(alpha_GC, p);
            vec6 = multi_scalaire(omega_GC, s);
            vec7 = Sommevect(vec5, vec6);
            x = Sommevect(x, vec7);

            vector<double> vec8(size);
            vec8 = multi_scalaire(omega_GC, t);
            r = Sousvect(s, vec8);

            norme_r = sqrt(produit_scalaire(r, r));
            norme_b = sqrt(produit_scalaire(b, b));

            iteration++;
        }

    return x;
}

int main(int argc, char **argv)
{
    // Données sur la boucle en temps
    int dt, nb_iterations;
    nb_iterations = 100;
    double tf = 10.0;
    dt = tf/nb_iterations;

    double alpha, beta, gamma;
    alpha = 1 + D * dt * (2. / pow(dx, 2) + 2. / pow(dy, 2));
    beta = -D * dt / pow(dx, 2);
    gamma = -D * dt / pow(dy, 2);

    double alpha_r, beta_r, coeff_Robin;
    alpha_r = 1. / 2.;
    beta_r = 1. / 2.;
    coeff_Robin = 2 * dt * beta_r * D / (alpha_r * dy);

    // Choix du cas test :
    int cas = 1;

    // Données pour le bi-gradient conjugué
    double erreur_schwartz = 1e-8;
    int iter = 0;
    int itermax = 10000;

    // Longueur du recouvrement 
    int taille_recouvrement = 1;


    vector<double> X(Nx);
    vector<double> stencil_up(3 * Nx), stencil_down(3 * Nx);

    for (int i=0 ; i<3*Nx; i++)
    {
        stencil_up[i] = 0.;
        stencil_down[i] = 0.;
    }

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int jbeg, jend;
    charge(me, Ny, np, &jbeg, &jend);

    MPI_Status Status;

    double debut = MPI_Wtime();

    for (int i = 0; i < Nx; i++)
    {
        X[i] = (i + 1) * dx;
    }

    int taille_0,taille_vect_0;
    vector<double> Y0,RHS0,U0;

    int taille_np,taille_vect_np;
    vector<double> Ynp,RHSnp,Unp;

    int taille_me,taille_vect_me;
    vector<double> Yme,RHSme,Ume;

    // Initialisation
    if (me == 0)
    {

        taille_0 = (jend - jbeg + 1) + taille_recouvrement;
        taille_vect_0 = taille_0 * Nx ;  
        Y0.resize(taille_0);
        RHS0.resize(taille_vect_0) ;
        U0.resize(taille_vect_0)   ;

        for (int j = 0; j < taille_0; j++)
        {

            Y0[j] = (j + 1) * dy;
          
        }

        for  (int k=0 ; k< taille_vect_0 ; k++){

            U0[k]= 0.;
        }
    }

    else if (me == np - 1)
    {

        taille_np = jend - jbeg + 1;
        taille_vect_np = taille_np*Nx;

        Ynp.resize(taille_np);
        Unp.resize(taille_vect_np);
        RHSnp.resize(taille_vect_np);

        for (int j = taille_np - 1; j >= 0; --j)
        {
            Ynp[j] = Ly - (j + 1) * dy;
                         
        }

        for (int k = 0 ; k<taille_vect_np ; k++){

            Unp[k]=0.; 

        }
    }

    else
    {

        taille_me = (jend - jbeg + 1) + taille_recouvrement + 1;
        taille_vect_me = taille_me*Nx;
        
        Yme.resize(taille_me);
        Ume.resize(taille_vect_me);
        RHSme.resize(taille_me);


        for (int j = 0; j < taille_me; j++)
        {
            Yme[j] = (jbeg + j) * dy;
            
        }

        for(int k=0; k<taille_vect_me; k++){
            Ume[k]=0.;

        }

    }

    double erreur = 1.;
    double tk;

    for(int k = 0; k < nb_iterations; k++)
    {
        tk = (k+1) * dt;
        while (iter<itermax && erreur < erreur_schwartz){
        
        if (me==0){

            
            //recv_stencil_down
            MPI_Recv(&stencil_down[0], 3*Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
            

            RHS0 = second_membre(beta,gamma,coeff_Robin,X,Y0,tk,stencil_up,stencil_down,me, taille_recouvrement,cas); //Y0 definie dans une boucle if (à resoudre) et ajout de la boucle du temps
            U0   = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHS0, me);
            for (int i=0;i<3*Nx;i++)
            {
                stencil_up[i]=U0[(taille_0-3)*Nx+i];   // taille_0 à definir ou recuperer
            }

            //send_U0 = send stencil_up
            MPI_Send(&stencil_up[0], 3*Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

        }

        else if (me==np-1)
        {
            
            //recv stencil_down
            MPI_Recv(&stencil_up[0], 3*Nx, MPI_DOUBLE, np-2, 0, MPI_COMM_WORLD, &Status);
            

            RHSnp = second_membre(beta,gamma,coeff_Robin,X,Ynp,tk,stencil_up,stencil_down,me, taille_recouvrement,cas); //Ynp definie dans une boucle if (à resoudre) et ajout de la boucle du temps
            Unp  = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHSnp, me);
            for(int i=0;i<3*Nx;i++)
            {
                stencil_down[i]=Unp[i];
            }

            //send U1
            MPI_Send(&stencil_down[0], 3*Nx, MPI_DOUBLE, np-2, 0, MPI_COMM_WORLD);

        }

        else{

            // recv stencil_up et  stencil_down
            MPI_Recv(&stencil_up[0], 3*Nx, MPI_DOUBLE, me+1, 0, MPI_COMM_WORLD, &Status);
            MPI_Recv(&stencil_down[0], 3*Nx, MPI_DOUBLE, me-1, 0, MPI_COMM_WORLD, &Status);
            
            
            RHSme = second_membre(beta,gamma,coeff_Robin,X,Ynp,tk,stencil_up,stencil_down,me, taille_recouvrement,cas); //idem
            Ume   = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHSnp, me);

            for (int i=0;i<3*Nx;i++){

                stencil_up[i]  = Ume [(taille_me-3)*Nx+i];
                stencil_down[i]= Ume [i];
                
            }

            // send Ume pour stencil_up et stencil_down
            MPI_Send(&stencil_down[0], 3*Nx, MPI_DOUBLE, me-1, 0, MPI_COMM_WORLD);
            MPI_Send(&stencil_down[0], 3*Nx, MPI_DOUBLE, me+1, 0, MPI_COMM_WORLD);

        }
        iter = iter + 1;
    }
    }
    
    double fin = MPI_Wtime();
    MPI_Finalize();

    return 0;
}
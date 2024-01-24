#include <vector>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
#include <fstream>

using namespace std;

int Nx = 10, Ny = 14;
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
        *ligne_beg = me * (m + 1) + 1;
        *ligne_end = (me + 1) * (m + 1) ;
    }
    else
    {
        *ligne_beg = r + me * m + 1;
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

//-------------------Fonction pour calculer l'erreur------------------------------------
  double normeL2(std::vector<double> vec1, std::vector<double> vec2)
  {

      double E = 0.;
      for (int i = 0; i < vec1.size(); ++i) {
          double diff = vec1[i] - vec2[i];
          E+= diff * diff;
      }

      return sqrt(E);
  }

//---------------------Produit matrice/vecteur----------------------------------

vector<double> Produit_matriciel_robin(int taille_recouvrement, double alpha, double beta, double gamma, double coeff_Robin, vector<double> U, int me)
{

    int ibeg, iend;
    charge(me, Ny, np, &ibeg, &iend);

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
        int Nme = (iend - ibeg + 1) + taille_recouvrement;
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

//---------------Fonction contenant les solutions exactes -------------------------------
  double sol_exa(double x, double y,int cas)
  {
    if (cas == 1) {
        return x*(1-x)*y*(1-y);
    }
    else if (cas == 2) {
        return sin(x)+cos(y);
    }

    else
    {
      return 0;
    }

  }

//---------------------vecteur second membre----------------------------------

vector<double> second_membre(double dt, double beta, double gamma, double coeff_Robin, vector<double> xcoord, vector<double> ycoord, double t, vector<double> stencil_up, vector<double> stencil_down, int me, vector<double> U, int taille_recouvrement, int cas)
{

    int sizex = xcoord.size();
    int sizey = ycoord.size();

    vector<double> S(sizex * sizey);

    for (int j = 1; j <= sizey; j++)
    {
        for (int i = 1; i <= sizex; i++)
        {
            S[(j - 1) * Nx + (i - 1)] = U[(j - 1) * Nx + (i - 1)] + dt * f(xcoord[i - 1], ycoord[j - 1], t, cas);
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
    int iterationmax = 100;
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

    r = Sousvect(b, Produit_matriciel_robin(taille_recouvrement, alpha, beta, gamma, coeff_Robin, x, me));
    rs = r;

    alpha_GC = 1.;
    omega_GC = 1.;

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
    int nb_iterations;
    double dt;
    nb_iterations = 100;
    double tf = 10.0;
    dt = tf/nb_iterations;
    //cout << "dt = " << dt << std::endl;

    double alpha, beta, gamma;
    alpha = 1 + D * dt * (2./pow(dx, 2) + 2./pow(dy, 2));
    beta = -D * dt/pow(dx, 2);
    gamma = -D * dt/pow(dy, 2);

    double alpha_r, beta_r, coeff_Robin;
    alpha_r = 1./2.;
    beta_r = 1./2.;
    coeff_Robin = 2 * dt * beta_r * D / (alpha_r * dy);

    // Choix du cas test :
    int cas = 3;

    // Données pour le bi-gradient conjugué
    double seuil_schwartz = 1e-8;
    int iter = 0;
    int itermax = 10;

    // Longueur du recouvrement 
    int taille_recouvrement = 1;


    vector<double> X(Nx);
    for (int i = 0; i < Nx; i++)
    {
        X[i] = (i + 1) * dx;
    }


//----------------------------------------Région MPI----------------------------------
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int jbeg, jend;
    charge(me, Ny, np, &jbeg, &jend);
    printf("me = %d , jbeg = %d, jend = %d\n", me, jbeg, jend);

    MPI_Status Status;

    double debut = MPI_Wtime();

    int taille_0,taille_vect_0;
    vector<double> Y0, RHS0, U0, prod_0;

    int taille_np,taille_vect_np;
    vector<double> Ynp, RHSnp, Unp, Unp_ex;

    int taille_me,taille_vect_me;
    vector<double> Yme, RHSme, Ume, Ume_ex;

    // Initialisation
    if (me == 0)
    {

        taille_0 = (jend - jbeg + 1) + taille_recouvrement;
        taille_vect_0 = taille_0 * Nx ;  
        Y0.resize(taille_0);
        RHS0.resize(taille_vect_0) ;
        U0.resize(taille_vect_0);
        prod_0.resize(taille_vect_0);

        for (int j = 0; j < taille_0; j++)
        {

            Y0[j] = (j + 1) * dy;
            //printf("me = %d, Y0[%d] = %f\n", me, j, Y0[j]);
          
        }

        for  (int k=0 ; k< taille_vect_0 ; k++){

            U0[k]= 1.;
        }

        prod_0 = Produit_matriciel_robin(taille_recouvrement, alpha, beta, gamma, coeff_Robin, U0, me);
        // for  (int k=0 ; k< taille_vect_0 ; k++){
        //     printf("me = %d, prod_0[%d] = %f\n", me, k, prod_0[k]);
        // }
    }

    else if (me == np - 1)
    {

        taille_np = jend - jbeg + 1;
        taille_vect_np = taille_np*Nx;

        Ynp.resize(taille_np);
        Unp.resize(taille_vect_np);
        Unp_ex.resize(taille_vect_np);
        RHSnp.resize(taille_vect_np);

        for (int j = taille_np - 1; j >= 0; --j)
        {
            Ynp[j] = Ly - (j + 1) * dy;
                         
        }

        for (int k = 0 ; k<taille_vect_np ; k++){

            Unp[k]=1.; 

        }
    }

    else
    {

        taille_me = (jend - jbeg + 1) + taille_recouvrement;
        taille_vect_me = taille_me*Nx;
        
        Yme.resize(taille_me);
        Ume.resize(taille_vect_me);
        Ume_ex.resize(taille_vect_me);
        RHSme.resize(taille_vect_me);


        for (int j = 0; j < taille_me; j++)
        {
            Yme[j] = (jbeg + j) * dy;
            //printf("me = %d, Yme[%d] = %f\n", me, j, Yme[j]);
            
        }

        for(int k=0; k<taille_vect_me; k++){
            Ume[k]=1.;

        }

    }

    vector<double> stencil_up(3 * Nx), stencil_down(3 * Nx);

    for (int i=0 ; i<3*Nx; i++)
    {
        stencil_up[i] = 1.;
        stencil_down[i] = 1.;
    }

    double erreur = 1.;
    double tk;

    double erreur_schwartz = 1.0;
    double err_scwhartz ;

    for(int k = 0; k < 1; k++)
    {
        tk = (k+1) * dt;
        //while (iter<itermax && erreur_schwartz > seuil_schwartz){
        while (iter<1 ){
        if (me==0){
            //double erreur_schwartz_0 = 0.;
    
            vector<double>Common_down((taille_recouvrement+1)*Nx);
         
            vector<double> U0_old(taille_vect_0);
            
            U0_old = U0 ;

            for (int i=0 ; i < (taille_recouvrement+1)*Nx ; i++)
            {
                Common_down[i] = U0_old[taille_vect_0-(taille_recouvrement+1)*Nx + i] ;
            }

            RHS0 = second_membre(dt, beta,gamma,coeff_Robin,X,Y0,tk,stencil_up,stencil_down,me, U0, taille_recouvrement,cas); 
            // for (int i=0 ; i<RHS0.size();i++)
            // {
            //     printf("me = %d, RHS0[%d] = %f\n", me, i, RHS0[i]);
            // }
            U0   = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHS0, me);
            for (int i=0 ; i<taille_vect_0;i++)
            {
                printf("U0[%d] = %f\n", i, U0[i]);
            }

            for (int i=0;i<3*Nx;i++)
            {
                stencil_up[i]=U0[(taille_0 - (taille_recouvrement + 2)) * Nx + i];   // taille_0 à definir ou recuperer
            }

            for (int i=0 ; i<(taille_recouvrement+1)*Nx ; i++)
            {
                err_scwhartz = fabs(Common_down[i] - U0[taille_vect_0-(taille_recouvrement+1)*Nx + i]);
            }

            MPI_Send(&stencil_up[0], 3*Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&stencil_up[0], 3*Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

        }

        else if (me==np-1)
        {
            //double erreur_schwartz_np = 0.;

            vector<double>Common_up_np((taille_recouvrement+1)*Nx);
         
            vector<double> Unp_old(taille_vect_np);
            
            Unp_old = Unp;

            for (int i=0 ; i < (taille_recouvrement+1)*Nx ; i++)
            {
                Common_up_np[i] = Unp_old[i] ;
            }

            RHSnp = second_membre(dt, beta,gamma,coeff_Robin,X,Ynp,tk,stencil_up,stencil_down,me, Unp, taille_recouvrement,cas); //Ynp definie dans une boucle if (à resoudre) et ajout de la boucle du temps
            Unp  = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHSnp, me);

            for(int i=0;i<3*Nx;i++)
            {
                stencil_down[i]=Unp[i + (taille_recouvrement - 1) * Nx];
            }

            for (int i=0 ; i<(taille_recouvrement+1)*Nx ; i++)
            {
                err_scwhartz = fabs(Common_up_np[i] - Unp[i]);
            }
            
            MPI_Send(&stencil_down[0], 3*Nx, MPI_DOUBLE, np-2, 0, MPI_COMM_WORLD);    
            MPI_Recv(&stencil_down[0], 3*Nx, MPI_DOUBLE, np-2, 0, MPI_COMM_WORLD, &Status);

        }

        else{
            //double erreur_schwartz_me = 0.;
            double err1,err2;

            vector<double>Common_down_me((taille_recouvrement+1)*Nx);
            vector<double>Common_up_me((taille_recouvrement+1)*Nx);

            vector<double> Ume_old(taille_vect_me);

            Ume_old = Ume ;

            for (int i=0 ; i < (taille_recouvrement+1)*Nx ; i++)
            {
                Common_down_me[i] = Ume[taille_vect_me-(taille_recouvrement+1)*Nx + i] ;
                Common_up_me[i]   = Ume[i] ;
            }

            RHSme = second_membre(dt, beta,gamma,coeff_Robin,X,Yme,tk,stencil_up,stencil_down,me, Ume, taille_recouvrement,cas); //idem
            Ume   = BiCGStab( taille_recouvrement, beta,alpha, gamma, coeff_Robin, RHSme, me);

            for (int i=0;i<3*Nx;i++){

                stencil_up[i]  = Ume [(taille_me - (taille_recouvrement + 2)) * Nx + i];
                stencil_down[i]= Ume [i + (taille_recouvrement - 1) * Nx];
                
            }

            for (int i=0 ; i<(taille_recouvrement+1)*Nx ; i++)
            {
                err1 = fabs(Ume[taille_vect_me-(taille_recouvrement+1)*Nx + i]-Common_down_me[i]);
                err2 = fabs(Ume[i]-Common_up_me[i]);

                err_scwhartz = fmax(err1, err2);
            }

            MPI_Send(&stencil_down[0], 3*Nx, MPI_DOUBLE, me-1, 0, MPI_COMM_WORLD);
            MPI_Send(&stencil_up[0], 3*Nx, MPI_DOUBLE, me+1, 0, MPI_COMM_WORLD);

            MPI_Recv(&stencil_up[0], 3*Nx, MPI_DOUBLE, me+1, 0, MPI_COMM_WORLD, &Status);
            MPI_Recv(&stencil_down[0], 3*Nx, MPI_DOUBLE, me-1, 0, MPI_COMM_WORLD, &Status);

        }
        MPI_Allreduce(&err_scwhartz,&erreur_schwartz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        iter = iter + 1;
    }

    if (me == 0)
    {   
        std::string filename = "Results_" + std::to_string(cas) + "_" + std::to_string(me) + "_" + std::to_string(k) +  ".dat";
        std::ofstream outfile(filename);
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < taille_0; j++)
            {
                outfile << X[i] << " " << Y0[j] << " " << U0[j * Nx + i ] << " " << sol_exa(X[i], Y0[j], cas) << std::endl;
            }
        }
        outfile.close();
    }

    else if (me == np - 1)
    {   
        std::string filename = "Results_" + std::to_string(cas) + "_" + std::to_string(me) + "_" + std::to_string(k) +  ".dat";
        std::ofstream outfile(filename);
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < taille_np; j++)
            {
                outfile << X[i] << " " << Ynp[j] << " " << Unp[j * Nx + i] << " " << sol_exa(X[i], Ynp[j], cas) << std::endl;
            }
        }
        outfile.close();
    }

    else 
    {
        std::string filename = "Results_" + std::to_string(cas) + "_" + std::to_string(me) + "_" + std::to_string(k) +  ".dat";
        std::ofstream outfile(filename);
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < taille_me; j++)
            {
                outfile << X[i] << " " << Yme[j] << " " << Ume[j * Nx + i ] << " " << sol_exa(X[i], Yme[j], cas) << std::endl;
            }
        }
        outfile.close();
    }

    }
    if (me == 0)
    {
        //Afficher U0
        //printf("U0:\n");
        for (int k = 0; k < taille_vect_0; k++) {
            //printf("U0[%d] = %f\n", k, U0[k]);
        }
        //printf("\n");

        // Afficher U0_ex
        //printf("U0_ex:\n");
        for (int k = 0; k < taille_vect_0; k++) { // Assurez-vous que taille_vect_0 est la bonne taille pour U0_ex
            //printf("U0_ex[%d] = %f\n", k, U0_ex[k]);
        }
        //printf("\n");

        double err_0;
        vector<double> U0_ex;
        U0_ex.resize(taille_vect_0);
        
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < taille_0; j++)
            {   
                int k;
                k = j * Nx + i;
                U0_ex[k] =  sol_exa(X[i], Y0[j], cas);
            }
        }
        err_0 = normeL2(U0,U0_ex);
        //printf("L'erreur est %.15f \n", err_0);
    }

    if (me == np - 1)
    {
        double err_np;
        vector<double> Unp_ex;
        Unp_ex.resize(taille_vect_np);
        
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < taille_np; j++)
            {   
                int k;
                k = j * Nx + i;
                Unp_ex[k] =  sol_exa(X[i], Ynp[j], cas);
            }
        }

        // Afficher Unp
        //printf("Unp:\n");
        for (int k = 0; k < taille_vect_np; k++) {
            //printf("Unp[%d] = %f\n", k, Unp[k]);
        }
        //printf("\n");

        // Afficher Unp_ex
        //printf("Unp_ex:\n");
        for (int k = 0; k < taille_vect_np; k++) { // Assurez-vous que taille_vect_0 est la bonne taille pour U0_ex
            //printf("Unp_ex[%d] = %f\n", k, Unp_ex[k]);
        }
        //printf("\n");

        err_np = normeL2(Unp,Unp_ex);
        //printf("L'erreur est %.15f \n", err_np);
    }
    
    double fin = MPI_Wtime();
    MPI_Finalize();

    return 0;
}
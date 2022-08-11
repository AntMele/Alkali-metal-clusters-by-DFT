#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int N = 40;
int n = 10000;
double h = 3e-3;
double * rho;
double * Vint;
double ** u;
double * mu;
double beta1 = 0.56371;
double beta2 = 0.27358;
double gamma = -0.103756;
double Rs = 3.93;
double Rc;
double alpha = 0.1;
int l = 0;
FILE * f1;
FILE * f2;

// Desmos graph at https://www.desmos.com/calculator/epri0j4xyf
double Vext(double r){
    if(r < Rc){
        return (pow(r,2) - 3*pow(Rc,2)) / (2*pow(Rs,3));
    }else{
        return -N / r;
    }
}

double epsilon_xc(double rho_i){
    double sqrt_rs = pow(3./(4*M_PI*rho_i), 1./6);
    return -3./4.*pow(3*rho_i/M_PI,1./3) + gamma/(1+beta1*sqrt_rs+beta2*pow(sqrt_rs,2));
}

double Vxc(double rho_i){
    double sqrt_rs = pow(3/(4*M_PI*rho_i), 1./6);
    return 2*pow(3*rho_i/M_PI,1./3) - gamma/3*(6+7*beta1*sqrt_rs+8*beta2*pow(sqrt_rs,2))/pow(1+beta1*sqrt_rs+beta2*pow(sqrt_rs,2),2);
}

double k_square(double E, int i){
    double r = h * i;
    return 2*E - l*(l+1)/pow(r,2) - 2*Vext(r) - 2*Vint[i] + Vxc(rho[i]);
}

double numerov(double E){
    double k_square0 = 0;
    double k_square1 = k_square(E, 1);
    double k_square2 = k_square(E, 2);
    double y0 = 0;
    double y1 = pow(h,l+1);
    double y2 = y1 * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 3; i <= n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, i);
        y0 = y1;
        y1 = y2;
        y2 = (y1 * (2-5*pow(h,2)*k_square1/6) - y0 * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }
    return y2;
}


double numerovSave(double E, double * y){
    double k_square0 = 0;
    double k_square1 = k_square(E, 1);
    double k_square2 = k_square(E, 2);
    y[0] = 0;
    y[1] = pow(h,l+1);
    y[2] = y[1] * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 2; i < n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, i+1);
        y[i+1] = (y[i] * (2-5*pow(h,2)*k_square1/6) - y[i-1] * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }

    // I correct the divergence of the wave function at rmax
    int index = 5000;
    double min = fabs(y[index]);
    for(int i = 5001; i <= n; i++){
        if(fabs(y[i]) < min){
            min = fabs(y[i]);
            index = i;
        }
    }
    for(int i = index+1; i <= n; i++){
        y[i] = min;
    }

    // I calculate the norm with the trapeziudal rule
    double norm = (pow(y[0],2) + pow(y[n],2)) / 2;
    for(int i = 1; i < n; i++){
        norm += pow(y[i], 2);
    }
    norm = sqrt(norm * h);
    for(int i = 0; i <= n; i++){
        y[i] /= norm;
    }

    return y[n];
}

double bisectionMethod(double xmin, double xmax, double f_xmin, double f_xmax){
    double f_x;
    double precision = 0.0000000000001;
    while(fabs((xmax - xmin)/(xmax + xmin)) >= precision){
        f_x = numerov((xmax+xmin)/2);
        if(f_x * f_xmin < 0){
            xmax = (xmax + xmin) / 2;
            f_xmax = f_x;
        }else{
            xmin = (xmax + xmin) / 2;
            f_xmin = f_x;
        }
    }

    return (xmax + xmin) / 2;
}


double zeroFinder(double xmin, double xmax, double dx){
    double f_xmin = numerov(xmin);
    double f_xmax;
    for(int i = 1; i <= (int)((xmax-xmin)/dx); i++){
        f_xmax = numerov(xmin + i*dx);
        if(f_xmin * f_xmax < 0){
            return bisectionMethod(xmin+(i-1)*dx, xmin+i*dx, f_xmin, f_xmax);
        }
        f_xmin = f_xmax;
    }
    printf("Error in zeroFinder: no zero found in the interval (%f, %f) with dx = %f\n", xmin, xmax, dx);
	exit(EXIT_FAILURE);
}

void computeVint(){
    double * vec1 = malloc((n+1)*sizeof(double));
    double * vec2 = malloc((n+1)*sizeof(double));
    vec1[0] = 0;
    vec2[n] = 0;
    for(int i = 1; i <= n; i++){
        vec1[i] = vec1[i-1] + rho[i]*pow(i,2);
        vec2[n-i] = vec2[n-i+1] + rho[n-i]*(n-i);
    }
    Vint[0] = 4*M_PI*pow(h,2)*vec2[0];
    for(int i = 1; i <= n; i++){
        Vint[i] = 4*M_PI*pow(h,2)*(vec1[i]/i + vec2[i]);
    }
    free(vec1);
    free(vec2);
}

double step(){
    // I find the minimum of the potential
    double Emin = Vext(0) + Vint[0] - 0.5*Vxc(rho[0]);
    fprintf(f1,"%.12f;", Emin);
    fprintf(f2,"%.12f;", rho[0]);
    for(int i = 1; i <= n; i++){
        double vi = Vext(i*h) + Vint[i] - 0.5*Vxc(rho[i]);
        if(vi < Emin){
            Emin = vi;
        }
        fprintf(f1,"%.12f;", vi);
        fprintf(f2,"%.12f;", rho[i]);
    }
    fprintf(f1,"\n"); fprintf(f2,"\n");

    l = 0;
    double * rho_old = malloc((n+1)*sizeof(double));
    switch(N){
    case 8:
        mu[0] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[0], u[0]);
        l = 1;
        mu[1] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[1], u[1]);

        // I update rho
        for(int i = 1; i <= n; i++){
            rho_old[i] = rho[i];
            rho[i] = alpha*(pow(u[0][i],2)+3*pow(u[1][i],2))/(2*M_PI*pow(i*h,2)) + (1-alpha)*rho[i];
        }
        rho_old[0] = rho[0];
        rho[0] = rho[1];
        break;
    case 20:
        mu[0] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[0], u[0]);
        mu[1] = zeroFinder(mu[0]*0.99, 0, fabs(mu[0]*0.99)/50);
        numerovSave(mu[1], u[1]);
        l = 1;
        mu[2] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[2], u[2]);
        l = 2;
        mu[3] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[3], u[3]);

        // I update rho
        for(int i = 1; i <= n; i++){
            rho_old[i] = rho[i];
            rho[i] = alpha*(pow(u[0][i],2)+pow(u[1][i],2)+3*pow(u[2][i],2)+5*pow(u[3][i],2))/(2*M_PI*pow(i*h,2)) + (1-alpha)*rho[i];
        }
        rho_old[0] = rho[0];
        rho[0] = rho[1];
        break;
    case 40:
        mu[0] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[0], u[0]);
        mu[1] = zeroFinder(mu[0]*0.99, 0, fabs(mu[0]*0.99)/50);
        numerovSave(mu[1], u[1]);
        l = 1;
        mu[2] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[2], u[2]);
        mu[3] = zeroFinder(mu[2]*0.99, 0, fabs(mu[2]*0.99)/50);
        numerovSave(mu[3], u[3]);
        l = 2;
        mu[4] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[4], u[4]);
        l = 3;
        mu[5] = zeroFinder(Emin, 0, fabs(Emin)/50);
        numerovSave(mu[5], u[5]);

        // I update rho
        for(int i = 1; i <= n; i++){
            rho_old[i] = rho[i];
            rho[i] = alpha*(pow(u[0][i],2)+pow(u[1][i],2)+3*(pow(u[2][i],2)+pow(u[3][i],2))+5*pow(u[4][i],2)+7*pow(u[5][i],2))/(2*M_PI*pow(i*h,2)) + (1-alpha)*rho[i];
        }
        rho_old[0] = rho[0];
        rho[0] = rho[1];
        break;
    }

    // I update Vint
    double * Vint_old = malloc((n+1)*sizeof(double));
    for(int i = 0; i <= n; i++){
        Vint_old[i] = Vint[i];
    }
    computeVint();

    // I compute gamma and the total energy
    double gammaj = 0;
    double Ej = 0;
    for(int i = 1; i <= n; i++){
        gammaj += pow(i,2) * rho[i] * (Vint[i] - Vint_old[i] - 0.5*(Vxc(rho[i]) - Vxc(rho_old[i])));
        Ej += pow(i,2) * rho[i] * (0.5*Vint[i] - Vint_old[i] + 0.5*Vxc(rho_old[i]) + epsilon_xc(rho[i]));
    }
    free(rho_old);
    free(Vint_old);
    gammaj *= 4*M_PI*pow(h,3);
    Ej *= 4*M_PI*pow(h,3);
    switch(N){
    case 8:
        Ej += 2 * (mu[0] + 3*mu[1]);
        break;
    case 20:
        Ej += 2 * (mu[0] + mu[1] + 3*mu[2] + 5*mu[3]);
        break;
    case 40:
        Ej += 2 * (mu[0] + mu[1] + 3*(mu[2] + mu[3]) + 5*mu[4] + 7*mu[5]);
        break;
    }
    // I compute the polarizability
    double deltaN = 0;
    for(int i = round(Rc/h); i <= n; i++){
        deltaN += pow(i,2) * rho[i];
    }
    printf("%.15f, %.15f, %.15f;\n", fabs(gammaj), Ej, pow(Rc,3)*(1+4*M_PI*pow(h,3)*deltaN/N));
    return gammaj;
}


// I assume the energy levels are ordered as 1s - 1p - 2s - 1d - 2p - 1f ...
int main(){
    Rc = Rs * pow(N,1./3);
    rho = malloc((n+1)*sizeof(double));
    Vint = malloc((n+1)*sizeof(double));
    f1 = fopen("f1.txt", "w"); // I write the potential on file
    f2 = fopen("f2.txt", "w"); // I write the density on file
    char str[50];
    FILE * file;
    switch(N){
    case 8:
        u = (double**)malloc(2*sizeof(double *));
        for(int i = 0; i < 2; i++){
            u[i] = (double*)malloc((n+1)*sizeof(double));
        }
        mu = malloc(2*sizeof(double));

        file = fopen("txt/Na/N8/E0_0.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[0][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N8/E0_1.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[1][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);

        for(int i = 1; i <= n; i++){
            rho[i] = alpha*(pow(u[0][i],2)+3*pow(u[1][i],2))/(2*M_PI*pow(i*h,2));
        }
        rho[0] = rho[1];
        break;
    case 20:
        u = (double**)malloc(4*sizeof(double *));
        for(int i = 0; i < 4; i++){
            u[i] = (double*)malloc((n+1)*sizeof(double));
        }
        mu = malloc(4*sizeof(double));

        file = fopen("txt/Na/N20/E0_0.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[0][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N20/E1_0.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[1][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N20/E0_1.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[2][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N20/E0_2.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[3][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);

        for(int i = 1; i <= n; i++){
            rho[i] = alpha*(pow(u[0][i],2)+pow(u[1][i],2)+3*pow(u[2][i],2)+5*pow(u[3][i],2))/(2*M_PI*pow(i*h,2));
        }
        rho[0] = rho[1];
        break;
    case 40:
        u = (double**)malloc(6*sizeof(double *));
        for(int i = 0; i < 6; i++){
            u[i] = (double*)malloc((n+1)*sizeof(double));
        }
        mu = malloc(6*sizeof(double));

        file = fopen("txt/Na/N40/E0_0.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[0][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N40/E1_0.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[1][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N40/E0_1.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[2][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N40/E1_1.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[3][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N40/E0_2.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[4][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);
        file = fopen("txt/Na/N40/E0_3.txt", "r");
        for(int i = 0; i <= n; i++){
            fgets(str, 50, file);
            u[5][i] = atof(strrchr(str, ';')+1);
        }
        fclose(file);

        for(int i = 1; i <= n; i++){
            rho[i] = alpha*(pow(u[0][i],2)+pow(u[1][i],2)+3*(pow(u[2][i],2)+pow(u[3][i],2))+5*pow(u[4][i],2)+7*pow(u[5][i],2))/(2*M_PI*pow(i*h,2));
        }
        rho[0] = rho[1];
        break;
    default:
        printf("are u stupid?\n");
        exit(EXIT_FAILURE);
        break;
    }

    /*for(int i = 0; i <= n; i++){
        rho[i] = N*exp(-pow(i*h/6,2))/pow(6*sqrt(M_PI),3);
    }*/

    computeVint();
    double gammaj;
    do{
        gammaj = step();
    }while(fabs(gammaj) > 0.000000000000001);

    for(int i = 0; i <= n; i++){
        fprintf(f1,"%.12f;", Vext(i*h) + Vint[i] - 0.5*Vxc(rho[i]));
        fprintf(f2,"%.12f;", rho[i]);
    }
    fprintf(f1,"\n");
    fprintf(f2,"\n");
    fclose(f1);
    fclose(f2);

    return 0;
}

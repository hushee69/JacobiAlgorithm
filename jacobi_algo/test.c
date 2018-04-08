#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>


#define MIN(a, b) (a<=b ? a : b)

// peut mettre des lignes et colonnes soit differentes soit meme
#define LIGNES					3
#define COLONNES				4

double **allocation_dynamic(int lignes, int colonnes);
double *allocation_dynamic_solution(int lig);
void get_matrix(double **matrix, int lig);
void get_b(double **matrix, int lig);
void affiche(double **matrix, int lig, int col);
double *jacobi_solutions(double **matrix, int lig, int no_of_iterations, double * sol, double eps);
void affiche_solution(double *matrix, int lig);
void free_matrix(double ***matrix, int lig);
void free_matrix_solution(double **matrix);
double calcul_erreur(double **matrix, double *solution, int lig);
void stabilite(double ** matrice, double delta, int lig);



 /* Matrice test */

void matrice_bord(double **matrix, int lig);
void matrice_franc(double **matrix, int lig);
void matrice_bord(double **matrix, int lig);
void matrice_ding_dong(double **matrix, int lig);
void matrice_hilbert(double **matrix, int lig);







double wiki_array[12]= // CE TABLEAU REMPLIS LA MATRICE PAR DES VALEAURS CHOISIS
{
	4, 1, 0, -1, 3, 6, -2, -5, -2, 1, 4, 6
};

double b_array[]=   // CE TABLEAU VA REMPLIR LA COLONNE B
{
	8, 3, 8, 5, 4, 8, 7, 10, 0, 1,3
};

double **allocation_dynamic(int lignes, int colonnes)
{
	// allocation de memoire pour les lignes
	double **a_retourner=(double**)malloc(sizeof(double*)*lignes);
	int i=0;

	for( i=0; i<lignes; ++i )
	{
		// allocation de memoire pour des colonnes
		a_retourner[i]=(double*)calloc(colonnes, sizeof(double));
	}

	return a_retourner;
}

double *allocation_dynamic_solution(int lig)
{
        	double *solution=(double*)calloc(lig, sizeof(double));

    return solution;
}

void get_matrix(double **matrix, int lig)
{
	int i=0, j=0, k=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=wiki_array[k];
			k++;
		}
	}

	return;
}

void get_b(double **matrix, int lig)
{
	int i=0, k=0;

	for( i=0; i<lig; ++i )
	{
		matrix[i][lig]=b_array[k];
		k++;
	}

	return;
}

void affiche(double **matrix, int lig, int col)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		printf("[");
		for( j=0; j<col; ++j )
		{
			printf(" %f ", matrix[i][j]);
		}

		printf("]\n");
	}
	printf("\n");

	return;
}

/* fonction qui calcule la norme euclidienne */
double norme2(double *x,int n){
	int i;
	double res=0;
	for(i=0;i<n;i++){
		res=res+x[i]*x[i];
	}
	return sqrt(res);
}


double *jacobi_solutions(double **matrix, int lig, int no_of_iterations, double * sol, double eps)
{
	int i=0, j=0,iters=0, a=0;
	int *k = &a;
	double *b=(double*)calloc(lig, sizeof(double));
    double *solution=(double*)calloc(lig, sizeof(double)); // j'ai fait un autre tableau qui va sauvgarde les ancien valeur de solution
    double *r=(double*)calloc(lig, sizeof(double));
    double x_y = 3.0;



    iters=(no_of_iterations<0 ? (no_of_iterations*-1) : no_of_iterations);// la valeur absolue de no_of_iterations
    double khra = 1;

	while( *k<iters && khra > eps)
	{
		for( i=0; i<lig; ++i )
		{
			b[i]=matrix[i][lig]; // remplis le vecteur b par les valeur de dernier colonne de la matrice cette colonne c'est le vecteur B
			for( j=0; j<lig; ++j )
			{
                if(i!=j){
					b[i]=(b[i]-(matrix[i][j]*sol[j]));}
            }

			solution[i] = b[i]/matrix[i][i];

            printf("solution %f\n", solution[i]);

            r[i] = (solution[i] - sol[i]);
		}
        printf("****\n");

		     // j'ai rajouté ça je suis pas suque c'est correct

        khra = norme2(r, lig);

        for(j=0; j<lig; j++)
            {
                sol[j] = solution[j];
            }
        *k++;

}
	if(*k == no_of_iterations)
	{
		printf("erreur de convergence de la methode de jacobi\n");
    }
	if( b) // c'est quoi ça??
	{
		free(b);
		b=NULL;
	}

	printf(" le nombre d'étiration utilisé est %d\n", iters);
	return solution;
}

void affiche_solution(double *matrix, int lig)
{
	int i=0;

	for( i=0; i<lig; ++i )
	{
		printf("[ x%d = %lf ]\n", (i+1), matrix[i]);
	}
	printf("\n");
}

void free_matrix(double ***matrix, int lig)
{
	int i=0;
	double **mat=*matrix;

	for( i=0; i<lig; ++i )
	{
		if( mat[i] )
		{
			free(mat[i]);
		}
	}

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}


void free_matrix_solution(double **matrix)
{
	double *mat=*matrix;

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}

double calcul_erreur(double **matrix, double *solution, int lig)
{
	double *ret_val=(double*)calloc(lig, sizeof(double)); // on a pas liberer la memoire ici !!!
	double *B=(double*)calloc(lig, sizeof(double));
	double erreur=0.0f;
	int i=0, j=0, k=0;

	for( i=0; i<lig; ++i )
	{
		B[i]=matrix[i][lig];
	}

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			ret_val[i]=0.0f;
			for( k=0; k<lig; ++k )
			{
				ret_val[i]=ret_val[i]+(matrix[i][k]*solution[k]);
				//printf("B[%d]=%f\n", k, B[k]);
				//printf("ret_val[%d]=%f\n", i, ret_val[i]);
			}
		}
	}

	for( i=0; i<lig; i++ )
	{
		ret_val[i]=ret_val[i]-B[i];
		erreur = erreur + pow(ret_val[i], 2);
	}
	erreur = sqrt(erreur);

	return erreur;
}

void stabilite(double ** matrice, double delta, int lig)
{
    int i=0;

	for( i=0; i<lig; ++i )
	{
		matrice[i][lig] = matrice[i][lig] + delta;
    }
}



        /*MATRICE TEST*/

void matrice_bord(double **matrix, int lig)
{
	int i=0;

	// je recopie ce code a partir de ton mail
	for( i=0; i<lig; ++i )
	{
		matrix[i][i]=1;
	}

	for( i=0; i<lig; ++i )
	{
		matrix[0][i]=pow(2.0f, (1-(i+1)));
		matrix[i][0]=matrix[0][i];
	}

	return;
}



void matrice_ding_dong(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=1/(2*lig-(i+1)-(j+1)+1.5);
		}
	}

	return;
}

void matrice_franc(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			if( (i+1)>=((j+1)+2) )
			{
				matrix[i][j]=0;
			}
			else
			{
				matrix[i][j]=MIN(i+1, j+1);
			}
		}
	}

	return;
}


void matrice_hilbert(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=(1.0f/((i+1)+(j+1)-1));
		}
	}

	return;
}

void matrice_kms(double **matrix, int lig, double p)
{
	if (p<= 0 || p>= 1)
	{
		printf("erreur");
	}
	else{
		int i=0, j=0;

		for( i=0; i<lig; ++i)
		{
			for( j=0; j<lig; ++j)
				{
				matrix[i][j] = pow( p ,fabs((i+1) - (j+1)) );
				}
		}
	}
}

void matrice_lehmer(double **matrix, int lig)
{
	int i=0, j=0;

		for( i=0; i<lig; ++i)
		{
			for( j=0; j<lig; ++j)
			{
				if (i<=j)
				{
					matrix[i][j] = (i+1/j+1);
				}
				else
				    {
					matrix[i][j] = (i/j);
				    }
}
}
}

void matrice_lotkin(double **matrix, int lig)
{
	int i=0, j=0;

		for( i=0; i<lig; ++i)
		{
			for( j=0; j<lig; ++j)
			{
				if( i == 1)
				{	matrix[i][j] = 1;
				}
				else
				{
					matrix[i][j] = (1.0f/((i+1) + (j+1) -1));
				}
}
}
}


void matrice_moler(double **matrix, int lig)
{
	int i=0, j=0;

		for( i=0; i<lig; ++i)
		{
			for( j=0; j<lig; ++j)
			{
				if( i == j)
				{
					matrix[i][j] = (double)(i+1);
				}
				else
				{
					matrix[i][j] = (MIN(i+1 ,j+1) - 2.0f);
				}
}
}
}


int main()
{
	double **mat=allocation_dynamic(LIGNES, COLONNES);
	double *solution=allocation_dynamic_solution(LIGNES);
	double eps = 0.01;
    double erreur_en_jacobi=0.0f;
    //double delta = 0.2; //pour la fonction stabilité
    clock_t debut=clock(); // ON COMMENCE à CONPTER LA VITESSE DE CONVERGENCE

	get_matrix(mat, LIGNES); // remplissage de la matrice par des valeurs choisis
    get_b(mat, LIGNES);// remplissage de la colonne B
	affiche(mat, LIGNES, COLONNES);// affichage de la matrice

	solution=jacobi_solutions(mat, LIGNES, 17, solution, eps);// le vecteur vide solution est remplis par les valeur solutions donné par la fonction jacobi_solutio
	affiche_solution(solution, LIGNES);// affichage des solutions

   /* stabilite(mat, delta,LIGNES);
    solution=jacobi_solutions(mat, LIGNES, 17, solution, eps);
    affiche(mat, LIGNES, COLONNES);
    affiche_solution(solution, LIGNES);*/


    free_matrix(&mat, LIGNES);                     /* libiré la matrice mat et reloué un espace mémoire pour elle */
    mat=allocation_dynamic(LIGNES, COLONNES);       /* et la remplire à nouveau pour utilisé la fonction calcul erreur */
    get_matrix(mat, LIGNES);
    get_b(mat, LIGNES);

    erreur_en_jacobi=calcul_erreur(mat, solution, LIGNES); // calculer l'eereur et l'affiché
    printf("erreur = %f\n", erreur_en_jacobi);

	free_matrix(&mat, LIGNES);          /* liberation de la mémoire */
	free_matrix_solution(&solution);

    clock_t fin=clock();
    double time_utilise=(double)(fin-debut)/CLOCKS_PER_SEC; // CALCULE DE TEMPS PASSER DANS L'EXECUTION

    printf("le temps utilise=%fs\n", time_utilise);

	getchar();
	return 0;
}
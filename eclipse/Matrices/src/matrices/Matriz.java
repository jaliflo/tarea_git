/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package matrices;

import java.awt.Dimension;
import java.util.Random;

/**
 *
 * @author galvez
 */
public class Matriz {
    private int[][]datos;
    private Random rnd = new Random();
    
    public Matriz(int filas, int columnas, boolean inicializarAleatorio){
        datos = new int[columnas][];
        for(int i=0; i<columnas; i++){
            datos[i] = new int[filas];
            if (inicializarAleatorio)
                for(int j=0; j<filas; j++)
                    datos[i][j] = rnd.nextInt(100);
        }
    }
    public Matriz(Dimension d, boolean inicializarAleatorio){
        this(d.height, d.width, inicializarAleatorio);
    }
    
    public Dimension getDimension(){
        return new Dimension(datos.length, datos[0].length);
    }
    
    public static Matriz sumarDosMatrices(Matriz a, Matriz b) throws DimensionesIncompatibles { 
        if(! a.getDimension().equals(b.getDimension())) throw new DimensionesIncompatibles("La suma de matrices requiere matrices de las mismas dimensiones");        
        int i, j, filasA, columnasA; 
        filasA = a.getDimension().height; 
        columnasA = a.getDimension().width; 
        Matriz matrizResultante = new Matriz(filasA, columnasA, false);
        for (j = 0; j < filasA; j++) { 
            for (i = 0; i < columnasA; i++) { 
                matrizResultante.datos[i][j] += a.datos[i][j] + b.datos[i][j]; 
            } 
        } 
        return matrizResultante; 
    } 
    
    public static Matriz invertirMatriz(Matriz a) throws DimensionesIncompatibles{
    	if(! (a.getDimension().height == a.getDimension().width)) throw new DimensionesIncompatibles("La suma de matrices requiere matrices de las mismas dimensiones");
        int dimension = a.getDimension().height;
        Matriz x = new Matriz(dimension, dimension, false);
        Matriz b = new Matriz(dimension, dimension, false);
        int indices[] = new int[dimension];
        for (int i=0; i<dimension; ++i) 
            b.datos[i][i] = 1;
 
 // Transform the matrix into an upper triangle
        gaus(a, indices);
 
 // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<dimension-1; ++i)
            for (int j=i+1; j<dimension; ++j)
                for (int k=0; k<dimension; ++k)
                    b.datos[indices[j]][k]
                    	    -= a.datos[indices[j]][i]*b.datos[indices[i]][k];
 
 // Perform backward substitutions
        for (int i=0; i<dimension; ++i) 
        {
            x.datos[dimension-1][i] = b.datos[indices[dimension-1]][i]/a.datos[indices[dimension-1]][dimension-1];
            for (int j=dimension-2; j>=0; --j) 
            {
                x.datos[j][i] = b.datos[indices[j]][i];
                for (int k=j+1; k<dimension; ++k) 
                {
                    x.datos[j][i] -= a.datos[indices[j]][k]*x.datos[k][i];
                }
                x.datos[j][i] /= a.datos[indices[j]][j];
            }
        }
        return x;
    }
 
// Method to carry out the partial-pivoting Gaussian
// elimination.  Here index[] stores pivoting order.
 
    public static void gaus(Matriz a, int indices[]) 
    {
        int n = indices.length;
        double c[] = new double[n];
 
 // Initialize the index
        for (int i=0; i<n; ++i) 
            indices[i] = i;
 
 // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i) 
        {
            double c1 = 0;
            for (int j=0; j<n; ++j) 
            {
                double c0 = Math.abs(a.datos[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }
 
 // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i) 
            {
                double pi0 = Math.abs(a.datos[indices[i]][j]);
                pi0 /= c[indices[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
 
   // Interchange rows according to the pivoting order
            int itmp = indices[j];
            indices[j] = indices[k];
            indices[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
                int pj = a.datos[indices[i]][j]/a.datos[indices[j]][j];
 
 // Record pivoting ratios below the diagonal
                a.datos[indices[i]][j] = pj;
 
 // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a.datos[indices[i]][l] -= pj*a.datos[indices[j]][l];
            }
        }
    }
    
    public static Matriz multiplicarDosMatrices(Matriz a, Matriz b) throws DimensionesIncompatibles {
    	if(!(a.getDimension().width == b.getDimension().height)) throw new DimensionesIncompatibles("Las columnas de A deben coincidir con las filas de B");
    	int filasA = a.getDimension().height;
    	int columnasA = b.getDimension().width;
    	Matriz matrizResultante = new Matriz(filasA, columnasA, false);
    	for(int i = 0; i < filasA; i++) {
    		for(int j = 0; j < columnasA; j++) {
    			int elemento = 0;
    			int filaDeB = 0;
    			for(int k = 0; k < columnasA; k++) {
    				elemento = elemento+(a.datos[k][i]*b.datos[filaDeB][j]);
    				filaDeB++;
    			}
    			matrizResultante.datos[i][j] = elemento;
    		}
    	}
    	
    	return matrizResultante;
    }

    @Override
    public String toString(){
        String ret = "";
        ret += "[\n";
        for (int i = 0; i < getDimension().width; i++) {
            ret += "(";
            for (int j = 0; j < getDimension().height; j++) {  
                ret += String.format("%3d", datos[i][j]); 
                if (j != getDimension().height - 1) ret += ", ";
            } 
            ret += ")";
            if (i != getDimension().width - 1) ret += ",";
            ret += "\n";
        } 
        ret += "]\n";
        return ret;
    }
}

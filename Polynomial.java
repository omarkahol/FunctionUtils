package com.okahol.FunctionUtils;
public class Polynomial extends Calculus
{
	
	
	public Vector polyfit(Vector x, Vector y, int degree, String method)
	{
		Matrix A = new Matrix(degree+1,"0");
		Vector b = new Vector(degree+1, "0");
	
		for(int i=0;i<degree+1;i++)
		{
			for(int j=0;j<degree+1;j++)
			{
				A.setElement(tracev(vpow(x, i+j)), i, j);
				
			}
			b.setElement(i, tracev(vprodv(vpow(x, i), y)));
			
		}

		Vector coeff = new Vector(degree+1,"0");
		if (method.equalsIgnoreCase("GAUSS"))
		{
			coeff = (solve(A, b,new Vector(degree+1,"0"),true));
			coeff.t();
		}
		else if (method.contains("SOR-"))
		{
			method.trim();
			double w = Double.parseDouble(method.substring(method.indexOf("-")+1, method.length()));
			coeff = (solvew(A, b,new Vector(degree+1,"0"),w,true));
			coeff.t();
		}
		else if (method.equalsIgnoreCase("INVERSE"))
		{
			A=inv(A,  false);
			coeff=vmprod(A, b);
		}
		else
		{
			System.out.println("Avaliable methods are 'SOR-omega_value' and 'GAUSS'. Please chosse one of them");
		}
		
		return coeff;	
	}
	
	public Vector vpolyval(Vector coeff, Vector xval)
	{
		
		Vector yval = xval.clone();
		
		for(int i=0;i<xval.length();i++)
		{
			double sum = 0;
			for(int k=0; k<coeff.length();k++)
			{
				sum += (coeff.getElements()[k]*Math.pow(xval.getElements()[i], coeff.length()-k-1));
			}
			yval.setElement(i, sum);
		}
		return yval;
	}
	
	double polyval(Vector coeff, double x0)
	{
		double sum = 0;
		for(int k=0; k<coeff.length();k++)
		{
			sum += (coeff.getElements()[k]*Math.pow(x0, coeff.length()-k-1));
		}
		return sum;
	}
}

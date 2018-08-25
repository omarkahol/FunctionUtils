package com.okahol.FunctionUtils;

public class Taylor 
{
	static int starter = 20;
	
	public static  double T_sin(double x)
	{
		double a = Math.floor(x/(2*Math.PI));
		a=a*2*Math.PI;
		if(starter==0)
		{
			starter=20;
			return x;
		}
		else
		{
			
			starter--;
			return (Math.pow(-1, starter+1)/FantasticFactorial.factorial(2*(starter+1)+1, 1))*Math.pow((x-a), 2*(starter+1)+1) + T_sin(x-a);
		}
	}
	public  static double T_cos(double x)
	{
		double a = Math.floor(x/(2*Math.PI));
		a=a*2*Math.PI;
		
		if(starter==0)
		{
			starter=20;
			return 1;
		}
		else
		{
			starter--;
			return (Math.pow(-1, starter+1)/FantasticFactorial.factorial(2*(starter+1), 1))*Math.pow((x-a), 2*(starter+1)) + T_cos(x-a);
		}
	}
}

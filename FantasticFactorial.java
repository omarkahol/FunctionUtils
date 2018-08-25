package com.okahol.FunctionUtils;
public class FantasticFactorial
{
public static double factorial(int n, int step)
{
	if(n<=step)
	{
		return n;
	}
	else
	{
		return n*factorial(n-step,step);
	}
}
public static void main(String[] args) {
	System.out.println(FantasticFactorial.factorial(20, 1));
}


}

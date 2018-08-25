package com.okahol.FunctionUtils;
import java.io.Serializable;
import java.util.Arrays;
public class Vector implements Serializable,Comparable<Object>
{
	private static final long serialVersionUID = 1L;
	private int length;
	private double[] elements;
	
	public Vector(int length, String type) 
	{
		this.length=length;
		elements=new double[this.length];
		if (type.equalsIgnoreCase("i"))
		{
			for (int i=0; i<this.length;i++)
			{
				this.elements[i]=i+1;
			}
		}
		else if (type.equalsIgnoreCase("d"))
		{
			for (int i=0; i<this.length;i++)
			{
				this.elements[i]=this.length-i;
			}
		}
		else if(type.contains("n:"))
		{
			try 
			{
				double value = Double.valueOf(type.substring(2).trim());
				for (int i=0; i<this.length;i++)
				{
					this.elements[i]=value;
				}
			}
			catch (Exception e)
			{
				for (int i=0; i<this.length;i++)
				{
					this.elements[i]=0;
				}
			}
		}
		else if(type.contains("l:"))
		{
			type=type.trim();
			type=type.substring(2);
			double start=Double.parseDouble(type.substring(0, type.indexOf(":",type.indexOf(":"))));
			type=type.substring(type.indexOf(":")+1,type.length());
			double end=Double.parseDouble(type);
			double step=(end-start)/this.length;
			for(int i=0; i<this.length;i++)
			{
				setElement(i, start+i*step);
			}
		}
	}
	public Vector(int length)
	{
		setLength(length);
	}
	
	public double[] getElements()
	{
		double[] theel=new double[this.length];
		for(int i=0; i<this.length;i++)
		{
			theel[i]=this.elements[i];
		}
		return theel;
	}
	public void setElements(double[] elements)
	{
		if(this.length==elements.length)
		{
			int count=0;
			for(double el : elements)
			{
				this.elements[count]=el;
				count++;
			}
		}
	}
	public void setLength(int length)
	{
		this.length=length;
		this.elements=new double[this.length];
	}
	public int length()
	{
		return this.length;
	}
	public void setElement(int position, double element)
	{
		this.elements[position]=element;
	}
	
	public void sum(Vector v)
	{
		if(this.length==v.length())
		{
			for(int i=0; i<this.length;i++)
			{
				setElement(i,this.elements[i]+v.getElements()[i]);
			}
		}
	}
	public void print(String format)
	{
		for(double k:this.elements)
		{
			System.out.print("|");
			System.out.printf(format,k);
			System.out.println("|");
		}
	}
	public void lambdaprod(double lambda)
	{
		for(int i=0; i<this.length;i++)
		{
			setElement(i,this.elements[i]*lambda);
		}
	}
	public void t()
	{
		for(int i=0;i<this.length;i++)
		{
			this.elements[i]=this.elements[this.length-i-1];
		}
	}
	
	public void pow(double degree)
	{
		for(int i =0; i<this.length;i++)
		{
			setElement(i,Math.pow(this.elements[i], degree));
		}
	}
	
	@Override 
	public Vector clone()
	{
		Vector v1 = new Vector(this.length);
		v1.setElements(this.getElements());
		return v1;
	}
	
	@Override
	public String toString()
	{
		return "Length: " + this.length() + "\nElements: "+this.elements.toString();
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof Vector && obj != null)
		{
			Vector v = (Vector) obj;
			
			if(v.length==this.length && Arrays.equals(this.elements,v.elements) )
			{
				return true;
			}
			return false;
		}
		return false;
	}
	@Override
	public int compareTo(Object o) {
		if (o != null && o instanceof Vector)
		{
			Vector v1 = (Vector) o;
			if (v1.length==this.length)
				return 0;
			else return (this.length>v1.length) ? 1 : -1;
		}
		else
			return -1;
	}
	


}

package com.okahol.FunctionUtils;
import java.io.Serializable;
public class Matrix implements Serializable,Comparable
{
	private static final long serialVersionUID = 1L;
	private int rows;
	private int col;
	private double[][] elements;
	
	
	
	
	public Matrix(int rows, int col, String type) 
	{	
		init(rows,col,type);
	}
	public Matrix (int rows, int col)
	{
		init(rows,col,"0");
	}
	public Matrix (int n, String type)
	{
		init(n,n,type);
	}
	
	
	
	
	private void init(int rows, int col, String type)
	{
		this.rows=rows;
		this.col=col;
		elements=new double[rows][col];
		
		if (type.equalsIgnoreCase("hilbert") && (rows==col))
		{
			for (int i=0; i<rows; i++)
			{
				for(int k = 0; k<rows;k++)
				{
					elements[i][k]=Math.pow((i+k+1),-1);
				}
				
			}
		}
		else if (type.contains("diag:")&&(rows==col))
		{
			String string="";
			for(int i=5;i<type.length();i++)
			{
				string+=Character.toString(type.charAt(i));
			}
			double el = Double.valueOf(string);
			for (int i=0; i<rows; i++)
			{
				elements[i][i]=el;
				
			}
		}
		else if (type.equalsIgnoreCase("eye")&&rows==col)
		{
			for (int i=0; i<rows; i++)
			{
				elements[i][i]=1;
				
			}
		}
		else 
		{
			try
			{
				double value = Double.valueOf(type);
				for (int i=0; i<rows; i++)
				{
					for(int k = 0; k<col;k++)
					{
						elements[i][k]=value;
					}
					
				}
			}
			catch (Exception e)
			{
				for (int i=0; i<rows; i++)
				{
					for(int k = 0; k<col;k++)
					{
						elements[i][k]=0;
					}
					
				}
			}
		}
	}
	
	public int getRow()
	{
		return this.rows;
	}
	
	public int getCol()
	{
		return this.col;
	}
	
	public double[][] getElements()
	{
		double[][] theel = new double[this.rows][this.col];
		for(int i=0; i<this.getRow();i++)
		{
			for(int k=0; k<this.getCol();k++)
			{
				theel[i][k]=this.elements[i][k];
			}
		}
		return theel;
	}
	public Vector getArow(int number)
	{
		Vector row = new Vector(this.col);
		for(int i=0;i<this.col;i++)
		{
			row.setElement(i, getElements()[number][i]);
		}
		return row;
	}
	public Vector getAcol(int number)
	{
		Vector col = new Vector(this.rows);
		for(int i=0;i<this.rows;i++)
		{
			col.setElement(i, getElements()[i][number]);
		}
		return col;
	}
	
	
	public void setElement(double element, int row, int col)
	{
		this.elements[row][col]=element;
	}
	public void setArow(Vector row, int number)
	{
		if(row.length()==this.col)
		{
			for(int i=0;i<this.col;i++)
			{
				setElement(row.getElements()[i],number,i);
			}
		}
	}
	public void setAcol(Vector col, int number)
	{
		if(col.length()==this.rows)
		{
			for(int i=0;i<this.rows;i++)
			{
				setElement(col.getElements()[i],i,number);
			}
		}
	}
	public void setAdiag(Vector diag, int number)
	{
		if(this.rows==this.col && (diag.length()+Math.abs(number))==this.rows)
		{
			for (int i=0; i<diag.length();i++)
			{
				if(number>=0)
				{
				this.elements[i][i+number]=diag.getElements()[i];
				}
				else
				{
					this.elements[i-number][i]=diag.getElements()[i];
				}
			}
		}
	}
	public void setElements(double[][] elements)
	{
		this.elements=elements;
	}
	
	public void print(String format)
	{
		for (int i=0; i<rows; i++)
		{
			for(int k = 0; k<col;k++)
			{
				System.out.print(" |");
				System.out.printf(format,elements[i][k]);
				System.out.print("| ");
			}
			System.out.println();
			System.out.println();
		}
	}
	public void sum(Matrix m)
	{
		if(m.rows==this.rows && m.col==this.col)
		{
			for(int i = 0; i<this.rows;i++)
			{
				for (int k=0;k<this.col;k++)
				{
					setElement(this.getElements()[i][k]+m.getElements()[i][k],i,k);
				}
			}
		}
	}
	public void lambda(double lambda)
	{
		for(int i=0;i<this.rows;i++)
		{
			for(int k=0;k<this.col;k++)
			{
				setElement(getElements()[i][k]*lambda,i,k);
			}
		}
	}
	@Override
	public String toString()
	{
		return "Rows->"+this.rows+" Columns->"+this.col;
	}
	
	@Override 
	public Matrix clone()
	{
		Matrix m = new Matrix(this.rows, this.col);
		m.setElements(this.getElements());
		return m;
	}
	@Override
	public int compareTo(Object obj) 
	{
		if (obj == null && !(obj instanceof Matrix))
			return -1;
		else {
			Matrix m = (Matrix) obj;
			int mprod= m.col*m.rows;
			int thisprod = this.col*this.rows;
			if(thisprod==mprod)
				return 0;
			else if(thisprod<mprod)
				return -1;
			else
				return 1;
			
		}
	}


}

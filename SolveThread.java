package com.okahol.FunctionUtils;
import java.util.concurrent.Callable;

public class SolveThread implements Callable<Vector>{
	private Matrix m;
	private Vector id;
	private Calculus c;
	public SolveThread(Matrix m, Vector id, Calculus c) {
		// TODO Auto-generated constructor stub
		this.m=m;
		this.id=id;
		this.c=c;
	}
	@Override
	public Vector call() {
		return c.solve(m, id, new Vector(id.length()), false);
	}

}

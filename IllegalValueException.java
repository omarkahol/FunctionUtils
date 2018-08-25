package com.okahol.FunctionUtils;

public class IllegalValueException extends RuntimeException
{
	public IllegalValueException()
	{
		super("The value chosen is illegal");
	}
	public IllegalValueException(String message)
	{
		super(message);
	}
}

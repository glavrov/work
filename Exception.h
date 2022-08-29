#pragma once
#ifndef _EXCEPTION_H
#define _EXCEPTION_H

#include <iostream>

class Exception
{
public:
	virtual void ShowMessage() = 0;
};

class BadSpaceDimensionException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nРазмерность пространства некорректна" << std::endl;
	}
};

class BadVectorDimensionException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nРазмерности векторов не совпадают" << std::endl;
	}
};

class ZeroDivideException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nДеление на 0" << std::endl;
	}
};

class LimitExceedingException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nПревышена размерность вектора" << std::endl;
	}
};

class PointCountOutOfRange : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nКоличество точек неверное" << std::endl;
	}
};

class BadSelectedMethodException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nНеверный выбор метода" << std::endl;
	}
};

class CircleRadiusOutOfRange : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\nРадиус окружности неверный" << std::endl;
	}
};
#endif // !_EXCEPTION_H
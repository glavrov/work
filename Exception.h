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
		std::cout << "\n����������� ������������ �����������" << std::endl;
	}
};

class BadVectorDimensionException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n����������� �������� �� ���������" << std::endl;
	}
};

class ZeroDivideException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n������� �� 0" << std::endl;
	}
};

class LimitExceedingException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n��������� ����������� �������" << std::endl;
	}
};

class PointCountOutOfRange : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n���������� ����� ��������" << std::endl;
	}
};

class BadSelectedMethodException : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n�������� ����� ������" << std::endl;
	}
};

class CircleRadiusOutOfRange : public Exception
{
public:
	void ShowMessage()
	{
		std::cout << "\n������ ���������� ��������" << std::endl;
	}
};
#endif // !_EXCEPTION_H
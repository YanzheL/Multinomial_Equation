//
//  main.cpp
//  Multinomial_Equation
//
//  Created by LI YANZHE on 08/03/2017.
//  Copyright © 2017 Yanzhe Lee. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>

#define UNIX


typedef std::vector<int> VArray_Int;
typedef std::vector<double> VArray_DB;

void Show_Index_Page()
{
    puts("\n--------------------------------------------------------------------------------");
    puts("|                                                                              |");
    puts("|         Copyright (C) 2016 Yanzhe Lee. All rights reserved.                  |");
    puts("|                                                                              |");
    puts("|                       Harbin Institute of Technology                         |");
    puts("|                                                                              |");
    puts("|         License GPLv3+: GNU GPL version 3 or later                           |");
    puts("|                                                                              |");
    puts("|         This is free software: you are free to change and redistribute it    |");
    puts("|                                                                              |");
    puts("|         Email: lee.yanzhe@yanzhe.org           Version 1.0.0                 |");
    puts("|                                                                              |");
    puts("--------------------------------------------------------------------------------");
    puts("|                        ---- Multinomial Equation ----                        |");
    puts("--------------------------------------------------------------------------------");
    puts("");
}

void Safe_Flush(FILE *fp)                                           //用于清空scanf缓冲区
{
    int ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n');
}

template <class VA>
void Remove_Same(VA &array,int k)
{
    VA arrayCopy;
    arrayCopy.assign(array.begin(),array.end());
    long i,j;
    double precision=pow(10, -k);
    for (i=0; i<arrayCopy.size(); ++i)
    {
        for (j=0; j<arrayCopy.size(); ++j)
        {
            if (i<arrayCopy.size()
                &&j<arrayCopy.size()
                &&fabs(arrayCopy.at(j)-arrayCopy.at(i))<=precision&&i!=j)
            {
                arrayCopy.erase(arrayCopy.begin()+j);
            }
        }
    }
    array.clear();
    array.assign(arrayCopy.begin(),arrayCopy.end());
    
}

class integer
{
private:
    int n=0;
public:
    bool integerFlag=false;
    
    integer(double xx)
    {
        if (Check_Integer(xx))
        {
            n=xx;
        }
    }
    
    bool Check_Integer(double x)                                    //检查传入的参数是否是整数
    {
        if (x==floor(x))
            integerFlag=true;
        return integerFlag;
    }
    
    VArray_Int factorize()                                          //因数分解
    {
        VArray_Int factors;
        factors.push_back(1);
        int i;
        int nCopy=abs(n);                                           //复制n，防止更改原始的n
        
        for (i=2; i<=nCopy; ++i)
        {
            if (nCopy%i==0)
            {
                int count=0;
                ++count;
                factors.push_back(i);                               //将不重复的因数放入factors中
                do
                {
                    nCopy=nCopy/i;
                    if (count>=2)
                    {
                        int temp=i;
                        temp*=i;                                    //如果有重复的因数则加入因数的重数次方
                        factors.push_back(temp);
                    }
                    ++count;
                }while (nCopy%i==0);
            }
        }
        if (n!=factors.back())
        {
            factors.push_back(n);
        }
        
        unsigned long beforeSize=factors.size();
        
        for (i=0; i<beforeSize; ++i)
        {
            int temp=factors.at(i);
            factors.push_back(-temp);                               //同时包含相应的负因数
        }
        
        return factors;                                             //返回存有所有因子的Vector
    }
    
};

class equation
{
private:
    int n=0;
    int i,j;
    VArray_DB coefficient;
    typedef struct posibleRoot
    {
        VArray_Int numerator;
        VArray_Int dominator;
    }posibleRoot;
public:
    bool integerFlag=true;
    equation(int nn)
    {
        n=nn;
    }
    void inputEfficient()
    {
        using namespace std;
        cout<<"Please input coefficients of the Multinomial Equation"<<endl;
        for (i=0; i<n+1; ++i)
        {
            double temp;
            while (scanf("%lf",&temp)!=1)
            {
                printf("\r");
                std::cout<<"Invalid input, please input again: ";
                Safe_Flush(stdin);
            }
            
            coefficient.push_back(temp);
            
            //            std::cout<<"Coefficient "<<i<<" = "<<coefficient.at(i)<<std::endl;
            
        }
        //        std::cout<<"Input END"<<std::endl;
    }
    double valueOfEqua(double x)                                    //代入x计算多项式的值
    {
        double sum=0;
        for (i=0; i<n+1; ++i)
        {
            sum+=coefficient.at(i)*pow(x, n-i);
        }
        return sum;
    }
    
    void checkIntegerCoeffi()                                       //检查是否是整系数多项式
    {
        for (i=0;i<n+1; i++)
        {
            integer n(coefficient.at(i));
            if(n.integerFlag==false)
            {
                integerFlag=false;
                break;
            }
        }
    }
    
    posibleRoot Struct_Possible_Root()                              //返回一个存有分子可能值和分母可能值的结构
    {
        posibleRoot posible;
        integer lowIndex(coefficient.back());
        integer highIndex(coefficient.front());
        posible.numerator=lowIndex.factorize();
        posible.dominator=highIndex.factorize();
        
        return posible;
    }
    
    VArray_DB Einstein()                                            //Einstein判别法
    {                                                               //如果一个整系数多项式有有理根，设此有理根为p/q，则
        posibleRoot posible=Struct_Possible_Root();                 //  p必然整除常数项系数
        VArray_DB finalResult;                                      //  q必然整除最高次项系数
        
        //        std::cout<<"Einstein Begin"<<std::endl;
        //        std::cout<<"posible.dominator.size() = "<<posible.dominator.size()<<std::endl;
        //        std::cout<<"posible.numerator.size() = "<<posible.numerator.size()<<std::endl;
        
        int idomi,jnum;
        for (idomi=0; idomi<posible.dominator.size(); ++idomi)
        {
            for (jnum=0; jnum<posible.numerator.size(); ++jnum)
            {
                double tempResult=(double)posible.numerator.at(jnum)/(double)posible.dominator.at(idomi);
                
                if (valueOfEqua(tempResult)==0)
                {
                    finalResult.push_back(tempResult);
                }
            }
        }
        
        return finalResult;
    }
    
    bool Root_Existance(double A,double B)
    {
        if (valueOfEqua(A)*valueOfEqua(B)<=0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    VArray_DB Bisection_Method(double intervalA,double intervalB,int accuracy)      //二分法，给定区间端点A、B
    {                                                                               //用于寻找无理根的近似解
        double middle=(intervalA+intervalB)/2;
        double precision=pow(10, -accuracy);
        int preciAccu=accuracy;
        
        if (n>=3&&accuracy>8)
        {
            preciAccu=8;
        }
        
        double equaPrecision=pow(10, -preciAccu);
        
        VArray_DB roots;
        if (fabs(valueOfEqua(intervalA))<=equaPrecision)
        {
            roots.push_back(intervalA);
            //            continueFlag=false;
            
        }
        else if (fabs(valueOfEqua(intervalB))<=equaPrecision)
        {
            roots.push_back(intervalB);
            //            continueFlag=false;8
        }
        
        
        //            printf("Value = %lf\n",fabs(valueOfEqua(middle)));
        while(true)
        {
            if ((intervalB-intervalA)<=precision/10)                                //如果区间长度小于精度
            {                                                                       //无论是否找到根都退出，防止无限循环
                break;
            }
            if (fabs(valueOfEqua(middle))<=equaPrecision)
            {
                roots.push_back(middle);                                            //把找到的根存入数组
                //                break;
            }
            
            
            if ((Root_Existance(intervalA, middle)==false)&&(Root_Existance(middle, intervalB)==false))
            {
                break;
            }
            
            if(Root_Existance(intervalA, middle))
            {
                //                if (n<=5)                                                           //如果3次方程以内，则查找所有根，否则只查找一个根
                //                {                                                                   //不然超过3次方程递归所需时间太长
                VArray_DB tempRootA=Bisection_Method(intervalA, middle, accuracy);
                long int sizeTA=tempRootA.size();
                if (sizeTA>0)
                {
                    Remove_Same(tempRootA, accuracy);
                    sizeTA=tempRootA.size();
                    for (i=0; i<sizeTA; ++i)
                    {
                        roots.push_back(tempRootA.at(i));
                    }
                    //                    break;
                }
                //                }
                //                else intervalB=middle;
                
            }
            
            if(Root_Existance(middle,intervalB))
            {
                //                if (n<=5)
                //                {
                VArray_DB tempRootB=Bisection_Method(middle, intervalB, accuracy);
                long int sizeTB=tempRootB.size();
                if (sizeTB>0)
                {
                    Remove_Same(tempRootB, accuracy);
                    sizeTB=tempRootB.size();
                    for (i=0; i<sizeTB; ++i)
                    {
                        roots.push_back(tempRootB.at(i));
                    }
                }
                //                }
                //                else intervalA=middle;
            }
            
            
            //            if (n<=5)
            //            {
            if (roots.size()>0)
                break;
            //            }
            //            else middle=(intervalA+intervalB)/2;
        }
        
        return roots;
    }
    
    
};



int main(int argc, const char * argv[])
{
    using namespace std;
#ifdef UNIX
    system("clear");
#endif
#ifdef WINDOWS
    system("cls");
#endif
    
    int n;
    int i;
    char preciseFlag;
    Show_Index_Page();
    cout<<"Please input maximum index of the Multinomial Equation : \n";
    while (scanf("%d",&n)!=1)
    {
        printf("\r");
        cout<<"Invalid input, please input again: ";
        Safe_Flush(stdin);
    }
    equation f(n);
    f.inputEfficient();
    f.checkIntegerCoeffi();
    if (f.integerFlag)
    {
        VArray_DB result=f.Einstein();
        if (result.size()>0)
        {
            cout<<"This equation exists rational root"<<endl;
            Remove_Same(result,3);
            unsigned long rootN=result.size();
            for (i=0; i<rootN; ++i)
            {
                cout<<result.at(i)<<"\t";
            }
            cout<<endl;
            if (rootN==n)
            {
                exit(0);
            }
        }
        else
        {
            cout<<"This equation do not exist rational root\n"<<endl;
        }
    }
    
    cout<<"Do you want to find the approximate root? Press y or n (Default y): ";
    Safe_Flush(stdin);
    scanf("%c",&preciseFlag);
    if (preciseFlag == '\n') preciseFlag = 'y';
    while ((preciseFlag != 'y'&&preciseFlag != 'n') || preciseFlag == '\n')
    {
        printf("\r");
        cout<<"Unavailable Choice, please choose again: ";
        Safe_Flush(stdin);
        if (scanf("%c", &preciseFlag) != 1)
        {
            puts("Input error");
            exit(1);
        }
    }
    cout<<endl;
    
    if (preciseFlag=='y')
    {
        double intervalA,intervalB;
        int accuracy;
        int displayAccu=accuracy;
        
        if (accuracy>7)
        {
            displayAccu=7;
        }
        
        cout<<"Please input the interval left endpoint A and right endpoint B"<<endl;
        while (scanf("%lf %lf",&intervalA,&intervalB)!=2)
        {
            printf("\r");
            cout<<"Invalid input, please input again: ";
            Safe_Flush(stdin);
        }
        if (intervalA>intervalB)
        {
            cout<<"Auto corrected the issue of A > B\n"<<endl;
            swap(intervalA,intervalB);
        }
        cout<<"Precision = 1x10^(-K), K is a integer no more than 10"<<endl;
        cout<<"Please input the precision index 'K': ";
        while (scanf("%d",&accuracy)!=1||accuracy>10)
        {
            printf("\r");
            cout<<"Invalid input, please input again: ";
            Safe_Flush(stdin);
        }
        cout<<endl;
        
        VArray_DB roots=f.Bisection_Method(intervalA, intervalB, accuracy);
        if (roots.size()>0)
        {
            cout<<"Congratulations! This equation exist roots as listed"<<endl;
            Remove_Same(roots,displayAccu);
            for (i=0; i<roots.size(); ++i)
            {
                printf("%.8lf\n",roots.at(i));
            }
        }
        else
            cout<<"This equation probaly do not exist root"<<endl;
    }
    
    return 0;
}

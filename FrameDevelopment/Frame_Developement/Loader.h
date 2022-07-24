#ifndef Loader_h
#define Loader_h

#include <iostream>
#include<Eigen/Dense>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<Eigen/StdVector>
#include <string.h>


class LoadFile
{
public:
    int m_rowNumber;
    int m_colNumber;
    std::string m_filename;
    LoadFile(std::string filename)
    {
        int rowNumber=0;
        int colNumber=0;
        std::fstream File;
        File.open(filename,std::ios::in);
        std::string line;
        if(File.is_open())
        {
            
            while(getline(File,line))
            {
                rowNumber++;
                std::stringstream ss;
                std::string str;
                ss<<line;
                getline(ss,str);
                for(int i=0;i<str.length();i++)
                {
                    if(str[i]==',')
                    {
                        colNumber++;
                    }
                }
            }
        }
        else
        {
            std::cout<<"Failed to Load File"<<'\n';
        }
        File.close();
        if(rowNumber==0)
        {
            colNumber=0;
        }
        else
        {
            colNumber=colNumber/rowNumber+1;
        }
        
        m_rowNumber=rowNumber;
        m_colNumber=colNumber;
        m_filename=filename;
        
        std::cout<<"Number of row= "<<m_rowNumber<<'\n';
        std::cout<<"Number of column= "<<m_colNumber<<'\n';
    }
    
    Eigen::MatrixXf CallMatrix()
    {
        Eigen::MatrixXf mat(m_rowNumber,m_colNumber);
        std::fstream myfile;
        myfile.open(m_filename,std::ios::in);
        if(myfile.is_open())
        {
            std::string mystring;
            int rowIndex=0;
            std::string line;
            while(getline(myfile,line))
            {
                std::stringstream ss(line);
                for(int i=0;i<m_colNumber;i++)
                {
                    getline(ss,mystring,',');
                    if(mystring=="Fixed")
                    {
                        mat(rowIndex,i)=0;
                    }
                    else if(mystring=="Free")
                    {
                        mat(rowIndex,i)=1;
                    }
                    else
                    {
                        mat(rowIndex,i)=std::stof(mystring);
                    }
                    
                }
                
                rowIndex++;
            }
        }
        else
        {
            std::cout<<"Failed to construct matrix"<<'\n';
        }
        myfile.close();
        
        return mat;
    }
    
   
};


#endif /* Loader_h */

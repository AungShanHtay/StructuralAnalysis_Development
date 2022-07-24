#ifndef EndMoment_h
#define EndMoment_h

#include <iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<Eigen/StdVector>
#include <string.h>
#include<algorithm>

template<typename Derived>
class EndMoment
{
public:
    Eigen::MatrixXf m_force;
    Eigen::MatrixXf m_element;
    Eigen::MatrixXf m_node;
    EndMoment (Eigen::MatrixBase<Derived> &force,Eigen::MatrixBase<Derived> &element,Eigen::MatrixBase<Derived> &node)
    {
        m_force=force;
        m_element=element;
        m_node=node;
    }
    
    Eigen::MatrixXf callEndReaction()
    {
        Eigen::MatrixXf Force(m_force.rows(),10);
        Force.col(6).setZero();
        Force.col(7).setZero();
        Force.col(8).setZero();
        Force.col(9).setZero();
        
        std::vector<float> element_col1(m_element.rows());
        std::vector<float> force_col2(m_force.rows());
        Eigen::VectorXf::Map(&element_col1[0], m_element.rows()) = m_element.col(0);
        Eigen::VectorXf::Map(&force_col2[0],m_force.rows())=m_force.col(1);
        
        for(int a=0;a<force_col2.size();a++)
        {
            
            auto it = std::find(element_col1.begin(), element_col1.end(), force_col2[a]);
            if (it != element_col1.end())
            {
                auto idx0=std::distance(element_col1.begin(),it);
                std::cout<<idx0<<"\n";
            }
        }
        
        for (int i = 0; i < m_force.rows(); i++)
        {
            
            Force(i,0)=i+1;
            Force(i,1)=m_force(i,1);
            // finding element number of force_col2 vector into m_element
            auto it = std::find(element_col1.begin(), element_col1.end(), force_col2[i]);
            
            if (it != element_col1.end())
            {
                auto idx=std::distance(element_col1.begin(),it);
                Force(i,2)=m_element(idx,1); //finding the index
                Force(i,3)=m_element(idx,2);
                Force(i,4)=m_force(i,3);
                Force(i,5)=m_force(i,4);
                
                
                
                if(m_force(i,2)==10)
                {
                    float P=0;
                    if(m_force(i,3)!=0 && m_force(i,4)!=0)
                    {
                        std::cout<<"Error 04: Point load cannot be applied twice"<<'\n';
                    }
                    if(m_force(i,3)!=0)
                    {
                        P=m_force(i,3);
                        
                    }
                    if(m_force(i,4)!=0)
                    {
                        P=m_force(i,4);
                    }
                    
                    float Li=m_force(i,5);
                    float Lj=m_force(i,6);
                    float L=Li+Lj;
                    float Ri=(P*pow(Lj,2)*(3*Li+Lj))/pow(L,3);
                    float Mi=(P*Li*pow(Lj,2))/pow(L,2);
                    float Rj=(P*pow(Li,2)*(Li+3*Lj))/pow(L,3);
                    float Mj=(P*Lj*pow(Li,2))/pow(L,2);
                    
                    Force(i,6)+=-Ri;
                    Force(i,7)+=-Mi;
                    Force(i,8)+=-Rj;
                    Force(i,9)+=-Mj;
                    
                }
                
                else if(m_force(i,2)==11)
                {
                    float w1=m_force(i,3);
                    float w2=m_force(i,4);
                    float Li=m_force(i,5);
                    float Lj=m_force(i,6);
                    float dX1=m_node((int)m_element(idx,1)-1,1);
                    float dX2=m_node((int)m_element(idx,2)-1,1);
                    float dY1=m_node((int)m_element(idx,1)-1,2);
                    float dY2=m_node((int)m_element(idx,2)-1,2);
                    float L=sqrt(pow(dX2-dX1,2)+pow(dY2-dY1,2));
                    
                    if(w1==w2)
                    {
                        float W=w1;
                        float a=Li;
                        float b=L-Lj;
                        float S1=a/L;
                        float S2=b/L;
                        float Mi=((W*pow(L,2)/12)*(6*pow(S2,2)-8*pow(S2,3)+3*pow(S2,4)))-((W*pow(L,2)/12)*(6*pow(S1,2)-8*pow(S1,3)+3*pow(S1,4)));
                        float Mj=((W*pow(L,2)/12)*(4*pow(S2,3)-3*pow(S2,4)))-((W*pow(L,2)/12)*(4*pow(S1,3)-3*pow(S1,4)));
                        float Ri=(Mi+((W*(b-a))*(0.5*(b-a)+(L-b)))-Mj)/L;
                        float Rj=(W*(b-a)-Ri);
                        Force(i,6)+=-Ri;
                        Force(i,7)+=-Mi;
                        Force(i,8)+=-Rj;
                        Force(i,9)+=Mj;
                    }
                    
                    else if(w1==0 && w2!=0)
                    {
                        float W=w2;
                        float a=Li;
                        float c=L-(Li+Lj);
                        float b=a+c;
                        float d=L-a-(2*c/3);
                        float A=pow(d,2)*(3*L-2*d);
                        float B=pow(c,2)/3;
                        float C=(0.5*L-b+(17*c/45));
                        float Ri=(W*c/(2*pow(L,3)))*(A-B*C);
                        float Rj=(0.5*W*c-Ri);
                        
                        float D=pow(d,2)*(d-L);
                        float E=(L/3)+(17*c/90)-(b/2);
                        float Mi=-(W*c/(2*pow(L,2)))*(D+B*E);
                        
                        float F=d*pow(d-L,2);
                        float G=pow(c,2)/6;
                        float H=(L/3)+(17*c/45)-b;
                        float Mj=(W*c/(2*pow(L,2)))*(F+G*H);
                        Force(i,6)+=-Ri;
                        Force(i,7)+=-Mi;
                        Force(i,8)+=-Rj;
                        Force(i,9)+=+Mj;
                    }
                    
                    else if(w2==0 && w1!=0)
                    {
                        float W=w1;
                        float a=Lj;
                        float c=L-(Li+Lj);
                        float b=a+c;
                        float d=L-a-(2*c/3);
                        float A=pow(d,2)*(3*L-2*d);
                        float B=pow(c,2)/3;
                        float C=(0.5*L-b+(17*c/45));
                        float Ri=(W*c/(2*pow(L,3)))*(A-B*C);
                        float Rj=(0.5*W*c-Ri);
                        
                        float D=pow(d,2)*(d-L);
                        float E=(L/3)+(17*c/90)-(b/2);
                        float Mi=-(W*c/(2*pow(L,2)))*(D+B*E);
                        
                        float F=d*pow(d-L,2);
                        float G=pow(c,2)/6;
                        float H=(L/3)+(17*c/45)-b;
                        float Mj=(W*c/(2*pow(L,2)))*(F+G*H);
                        Force(i,6)+=-Rj;
                        Force(i,7)+=-Mj;
                        Force(i,8)+=-Ri;
                        Force(i,9)+=+Mi;
                    }
                    
                    else if(abs(w1)<abs(w2) && w1!=0 && w2!=0)
                    {
                        float W1=w1;
                        float a1=Li;
                        float b1=L-Lj;
                        float S1=a1/L;
                        float S2=b1/L;
                        float Mi1=((W1*pow(L,2)/12)*(6*pow(S2,2)-8*pow(S2,3)+3*pow(S2,4)))-((W1*pow(L,2)/12)*(6*pow(S1,2)-8*pow(S1,3)+3*pow(S1,4)));
                        float Mj1=((W1*pow(L,2)/12)*(4*pow(S2,3)-3*pow(S2,4)))-((W1*pow(L,2)/12)*(4*pow(S1,3)-3*pow(S1,4)));
                        float Ri1=(Mi1+((W1*(b1-a1))*(0.5*(b1-a1)+(L-b1)))-Mj1)/L;
                        float Rj1=(W1*(b1-a1)-Ri1);
                        
                        float W2=(w2-w1);
                        float a2=Li;
                        float c2=L-(Li+Lj);
                        float b2=a2+c2;
                        float d2=L-a2-(2*c2/3);
                        float A=pow(d2,2)*(3*L-2*d2);
                        float B=pow(c2,2)/3;
                        float C=(0.5*L-b2+(17*c2/45));
                        float Ri2=(W2*c2/(2*pow(L,3)))*(A-B*C);
                        float Rj2=(0.5*W2*c2-Ri2);
                        float D=pow(d2,2)*(d2-L);
                        float E=(L/3)+(17*c2/90)-(b2/2);
                        float Mi2=-(W2*c2/(2*pow(L,2)))*(D+B*E);
                        float F=d2*pow(d2-L,2);
                        float G=pow(c2,2)/6;
                        float H=(L/3)+(17*c2/45)-b2;
                        float Mj2=(W2*c2/(2*pow(L,2)))*(F+G*H);
                        
                        float Ri3=Ri1+Ri2;
                        float Mi3=Mi1+Mi2;
                        float Rj3=Rj1+Rj2;
                        float Mj3=Mj1+Mj2;
                        
                        Force(i,6)+=-Ri3;
                        Force(i,7)+=-Mi3;
                        Force(i,8)+=-Rj3;
                        Force(i,9)+=+Mj3;
                        
                    }
                    
                    else if(abs(w1)>abs(w2) && w1!=0 && w2!=0)
                    {
                        float W1=w2;
                        float a1=Li;
                        float b1=L-Lj;
                        float S1=a1/L;
                        float S2=b1/L;
                        float Mi1=((W1*pow(L,2)/12)*(6*pow(S2,2)-8*pow(S2,3)+3*pow(S2,4)))-((W1*pow(L,2)/12)*(6*pow(S1,2)-8*pow(S1,3)+3*pow(S1,4)));
                        float Mj1=((W1*pow(L,2)/12)*(4*pow(S2,3)-3*pow(S2,4)))-((W1*pow(L,2)/12)*(4*pow(S1,3)-3*pow(S1,4)));
                        float Ri1=(Mi1+((W1*(b1-a1))*(0.5*(b1-a1)+(L-b1)))-Mj1)/L;
                        float Rj1=(W1*(b1-a1)-Ri1);

                        float W2=(w1-w2);
                        float a2=Lj;
                        float c2=L-(Li+Lj);
                        float b2=a2+c2;
                        float d2=L-a2-(2*c2/3);
                        float A=pow(d2,2)*(3*L-2*d2);
                        float B=pow(c2,2)/3;
                        float C=(0.5*L-b2+(17*c2/45));
                        float Ri2=(W2*c2/(2*pow(L,3)))*(A-B*C);
                        float Rj2=(0.5*W2*c2-Ri2);
                        float D=pow(d2,2)*(d2-L);
                        float E=(L/3)+(17*c2/90)-(b2/2);
                        float Mi2=-(W2*c2/(2*pow(L,2)))*(D+B*E);
                        float F=d2*pow(d2-L,2);
                        float G=pow(c2,2)/6;
                        float H=(L/3)+(17*c2/45)-b2;
                        float Mj2=(W2*c2/(2*pow(L,2)))*(F+G*H);
                        
                        float Ri3=Ri1+Rj2;
                        float Mi3=Mi1+Mj2;
                        float Rj3=Rj1+Ri2;
                        float Mj3=Mj1+Mi2;
                        
                        Force(i,6)+=-Ri3;
                        Force(i,7)+=-Mi3;
                        Force(i,8)+=-Rj3;
                        Force(i,9)+=+Mj3;
                        
                    }
                    
                }
                
                
            }
            else
            {
                std::cout << "element member does not match with element table" << '\n';
            }
        }
        
        return Force;
    }
    
    
};



#endif /* EndMoment_h */

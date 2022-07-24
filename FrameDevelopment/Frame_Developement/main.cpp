
#include <iostream>
#include<Eigen/Dense>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<vector>
#include<iomanip>
#include<Eigen/StdVector>
#include "Loader.h"
#include "matplotlibcpp.h"
#include "EndMoment.h"

namespace plt=matplotlibcpp;

int main()
{
    
    LoadFile nodefile("Node.txt");
    Eigen::MatrixXf node=nodefile.CallMatrix();
    
    //**************Node table************************************************************************
    int Nwidth=14; //node table width
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Node table"<< "\n";
    std::cout << "No" << std::setw(Nwidth) << "X-Value" << std::setw(Nwidth)<<"Y-Value"<< '\n';
    for (int i = 0; i < node.rows(); i++)
    {
        std::cout << node(i, 0) << std::setw(Nwidth) << node(i, 1) << std::setw(Nwidth) << node(i, 2) << '\n';
    }
    
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";
    
    //**************Element table************************************************************************
    
    LoadFile elementfile("Element.txt");
    Eigen::MatrixXf element=elementfile.CallMatrix();
    int Ewidth=20;
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Element table"<< "\n";
    std::cout << "No" << std::setw(Ewidth) << "ithNode(m)" << std::setw(Ewidth)<<"jthNode(m)"<< std::setw(Ewidth)<<"E(kN/m2)"<<std::setw(Ewidth)<<"I(m4)"<<'\n';
    
    for (int i = 0; i < element.rows(); i++)
    {
        std::cout << element(i, 0) << std::setw(Ewidth) << element(i, 1) << std::setw(Ewidth) << element(i, 2) <<std::setw(Ewidth) << element(i, 3)<< std::setw(Ewidth) << element(i, 4)<< '\n';
    }
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";
    
    //**************Stiffness Matrix*************************************************************
    int dof=2*(int)node.rows();
    Eigen::MatrixXf K=Eigen::MatrixXf::Zero(dof,dof);
    
    for(int i=0;i<element.rows();i++)
    {
        int posA=(int)element(i,1)*3-(int)element(i,1)-2;
        int posB=(int)element(i,1)*3-(int)element(i,1)-1;
        int posC=(int)element(i,2)*3-(int)element(i,2)-2;
        int posD=(int)element(i,2)*3-(int)element(i,2)-1;
        float deltaX=pow(abs(node((int)element(i,1)-1,1)-node((int)element(i,2)-1,1)),2);
        float deltaY=pow(abs(node((int)element(i,1)-1,2)-node((int)element(i,2)-1,2)),2);
        float L=sqrt(deltaX+deltaY);
        float E=element(i,3);
        float I=element(i,4);
        
        std::cout<<"Element "<<i+1<<'\n';
        std::cout<<"Node"<<element(i,1)<<"("<<posA<<","<<posB<<")"<<"<---connected to------>"<<"Node"<<element(i,2)<<"("<<posC<<","<<posD<<")"<<'\n';
        std::cout<<"Length="<<L<<"(m)"<<'\n';
        std::cout<<"Modulus of Elasticity="<<E<<"(kN/m2)"<<'\n';
        std::cout<<"Moment of Inertia="<<I<<"(m4)"<<"\n\n";
        
        K(posA,posA)+=12*E*I/pow(L,3);
        K(posA,posB)+=6*E*I/pow(L,2);
        K(posA,posC)+=-12*E*I/pow(L,3);
        K(posA,posD)+=6*E*I/pow(L,2);
        
        K(posB,posA)+=6*E*I/pow(L,2);
        K(posB,posB)+=4*E*I/pow(L,1);
        K(posB,posC)+=-6*E*I/pow(L,2);
        K(posB,posD)+=2*E*I/pow(L,1);
        
        K(posC,posA)+=-12*E*I/pow(L,3);
        K(posC,posB)+=-6*E*I/pow(L,2);
        K(posC,posC)+=12*E*I/pow(L,3);
        K(posC,posD)+=-6*E*I/pow(L,2);
        
        K(posD,posA)+=6*E*I/pow(L,2);
        K(posD,posB)+=2*E*I/pow(L,1);
        K(posD,posC)+=-6*E*I/pow(L,2);
        K(posD,posD)+=4*E*I/pow(L,1);
        
    }
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Global Stiffness Matrix K"<< "\n";
    std::cout<<K<<"\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";
    
    //**************Support Conditions************************************************************************
    LoadFile supportconditionfile("SupportConditions.txt");
    Eigen::MatrixXf support=supportconditionfile.CallMatrix();
    int supportwidth=16;
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Support Conditions"<< "\n";
    std::cout << "Node" << std::setw(supportwidth) << "dY(m)" << std::setw(supportwidth)<<"θz(Degree)"<<'\n';
    
    for (int i = 0; i < support.rows(); i++)
    {
        std::cout << support(i, 0) << std::setw(supportwidth) << support(i, 1) << std::setw(supportwidth) << support(i, 2)<< '\n';
    }
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";
    //**************Load Conditions & Calculation of End Moment*********************************************************
    LoadFile forcefile("ForceConditions.txt");
    Eigen::MatrixXf force=forcefile.CallMatrix();
    int Fwidth=13; //node table width    
    EndMoment endmoment(force,element,node);
    Eigen::MatrixXf Force=endmoment.callEndReaction();
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Force table and Fixed End Reactions"<< "\n";
    std::cout << "No" << std::setw(Fwidth) << "Element" << std::setw(Fwidth)<<"Node ith"<< std::setw(Fwidth)<<"Node jth"<<std::setw(Fwidth)<<"Load1"<<std::setw(Fwidth)<<std::setw(Fwidth)<<"Load2"<<std::setw(Fwidth)<<"Ri"<<std::setw(Fwidth)<<"Mi"<<std::setw(Fwidth)<<"Rj"<<std::setw(Fwidth)<<"Mj"<<'\n';
    for (int i = 0; i < force.rows(); i++)
    {
        std::cout << Force(i, 0) << std::setw(Fwidth) << Force(i, 1) << std::setw(Fwidth) << Force(i, 2) <<std::setw(Fwidth) << Force(i, 3) <<std::setw(Fwidth) << Force(i, 4) <<std::setw(Fwidth) << Force(i, 5) <<std::setw(Fwidth) << Force(i, 6)  <<std::setw(Fwidth) << Force(i, 7) <<std::setw(Fwidth) << Force(i, 8)  <<std::setw(Fwidth) << Force(i, 9)<<'\n';
    }
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";

    //**************Boundary Condition************************************************************************
    
    Eigen::MatrixXf Boundary(node.rows(),5);
    Boundary.col(0)=node.col(0);
    Boundary.col(1)=support.col(1);
    Boundary.col(2)=support.col(2);
    Boundary.col(3).setZero();
    Boundary.col(4).setZero();
    int BDwidth=16;
    
    for(int i=0;i<Boundary.rows();i++)
    {
        for(int j=0;j<Force.rows();j++)
        {
            if(Boundary(i,0)==Force(j,2))
            {
                Boundary(i,3)+=-Force(j,6);
                Boundary(i,4)+=-Force(j,7);
            }
            
            if(Boundary(i,0)==Force(j,3))
            {
                Boundary(i,3)+=-Force(j,8);
                Boundary(i,4)+=-Force(j,9);
            }
            
        }
    }
        
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "Boundary Conditions"<< "\n";
    std::cout << "Node" << std::setw(BDwidth) << "dY(m)" << std::setw(BDwidth)<<"θz(Degree)"<<std::setw(BDwidth)<<"Fy(kN)"<<std::setw(BDwidth)<<"Mz(kN-m)"<<'\n';
    
    for (int i = 0; i < Boundary.rows(); i++)
    {
        std::cout << Boundary(i, 0) << std::setw(BDwidth) << Boundary(i, 1) << std::setw(BDwidth) << Boundary(i, 2) <<std::setw(BDwidth) << Boundary(i, 3) <<std::setw(BDwidth) << Boundary(i, 4) << '\n';
    }
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"<< "\n";
    std::cout << "\n\n";
    
    //**************Displacement & ForceVector********************************************************
    
    //Equivalent Nodal Forces Vector Fe
    Eigen::VectorXf Fequi(2*Boundary.rows());
    int countFe=0;
    for(int i=0;i<Boundary.rows();i++)
    {
        Fequi(countFe)=Boundary(i,3);
        countFe++;
        Fequi(countFe)=Boundary(i,4);
        countFe++;
    }
    
    std::cout<<"Equivalent Nodal Forces"<<'\n';
    std::cout<<Fequi<<"\n\n";
    
    Eigen::VectorXf Utemp=Eigen::VectorXf::Ones(dof);
    Eigen::VectorXf Ftemp=Eigen::VectorXf::Zero(dof);
    int count1=0;
    int count2=0;
    int KnownDisplacement=0; //=Unknown force
    
    for(int i=0;i<Boundary.rows();i++)
    {
        Utemp(count1)=Boundary(i,1);
        count1++;
        Utemp(count1)=Boundary(i,2);
        count1++;
        Ftemp(count2)=Boundary(i,3);
        count2++;
        Ftemp(count2)=Boundary(i,4);
        count2++;
    }
    
    for(int i=0;i<Utemp.rows();i++)
    {
        if(Utemp(i)!=1)
        {
            KnownDisplacement++;
        }
    }
    
    
    int KnownForce=dof-KnownDisplacement; //=Unknown Displacement
    //Known displacement and displacement vector
    Eigen::VectorXf UE(KnownDisplacement);
    Eigen::VectorXf FF(KnownForce);
    count1=0;
    count2=0;
    for(int i=0;i<Utemp.rows();i++)
    {
        if(Utemp(i)==0)
        {
            UE(count1)=Utemp(i);
            count1++;
        }
        else
        {
            FF(count2)=Ftemp(i);
            count2++;
        }
        
    }
    
    int unknownDisplacement=dof-KnownDisplacement;
    int unknownForce=dof-KnownForce;
    
    Eigen::MatrixXf KE(unknownForce, KnownDisplacement);
    Eigen::MatrixXf KEF(unknownForce, unknownDisplacement);
    Eigen::MatrixXf KF(KnownForce, unknownDisplacement);
    Eigen::MatrixXf KFE(unknownDisplacement, unknownForce);
    
    int rowindex1 = 0;
    int rowndex2 = 0;
    for (int i = 0; i < Utemp.rows(); i++)
    {
        int colindex1 = 0;
        int colindex2 = 0;
        if (Utemp(i) == 0)
        {
            // UE(i)=Utemp(i);
            for (int j = 0; j < Utemp.rows(); j++)
            {
                
                if (Utemp(j) == 0)
                {
                    KE(rowindex1, colindex1) = K(i, j);
                    colindex1++;
                }
                else
                {
                    KEF(rowindex1, colindex2) = K(i, j);
                    colindex2++;
                }
            }
            rowindex1++;
        }
        else
        {
            for (int k = 0; k < Utemp.rows(); k++)
            {
                if (Utemp(k) != 0)
                {
                    KF(rowndex2, colindex2) = K(i, k);
                    colindex2++;
                }
            }
            rowndex2++;
        }
    }
    
    KFE = KEF.transpose();
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Known Displacement Matrix KE" << '\n';
    std::cout << KE << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Matrix KEF" << '\n';
    std::cout << KEF << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Matrix KFE" << '\n';
    std::cout << KFE << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Matrix KF" << '\n';
    std::cout << KF << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    //--------------Displacement vector & Force vector-----------
    
    Eigen::VectorXf UF(KnownForce);
    UF = KF.colPivHouseholderQr().solve(FF);
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Solutions for unknown displacement vectors "
    << "\n";
    std::cout << UF << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    Eigen::VectorXf FE(KnownDisplacement);
    FE = KEF * UF;
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Solution for unknown force vectors "
    << "\n";
    std::cout << FE << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    Eigen::VectorXf Unew(dof);
    Eigen::VectorXf Feff(dof);
    Eigen::VectorXf Fnew(dof);
    // Unew<<UE,Uf;
    int indexE = 0;
    int indexf = 0;
    for (int i = 0; i < dof; i++)
    {
        if (Utemp(i) == 0)
        {
            Unew(i) = UE(indexE);
            Feff(i) = FE(indexE);
            indexE++;
        }
        else
        {
            Unew(i) = UF(indexf);
            Feff(i) = FF(indexf);
            indexf++;
        }
    }
    
    Fnew=Feff-Fequi;
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "New Displacement Vector"
    << "\n";
    std::cout << Unew << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "Effective Force Vector"
    << "\n";
    std::cout << Feff << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "New Force Vector"
    << "\n";
    std::cout << Fnew << "\n\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    
    Eigen::MatrixXf Knew(dof, dof);
    Knew.block(0, 0, KnownDisplacement, KnownDisplacement) = KE;
    Knew.block(0, KnownDisplacement, KnownDisplacement, KnownForce) = KEF;
    Knew.block(KnownDisplacement, 0, KnownForce, KnownDisplacement) = KFE;
    Knew.block(KnownDisplacement, KnownDisplacement, KnownForce, KnownForce) = KF;
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n";
    std::cout << "New Stiffness Matrix"
    << "\n";
    std::cout << Knew << "\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------------------"
    << "\n\n";
    
}

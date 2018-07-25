#include "MieosProjectOrganiserExtLib.hpp"
#include <iostream>

MieosProjectOrganiserExtLib::MieosProjectOrganiserExtLib(){

}

MieosProjectOrganiserExtLib::~MieosProjectOrganiserExtLib(){

}

void MieosProjectOrganiserExtLib::staticSayHi(){
   std::cout <<"Static Hi from external example class"<<std::endl;
}

void MieosProjectOrganiserExtLib::sayHi(){ 
   std::cout <<"Hi from external example class"<<std::endl;
}

#include <iostream>
#include <fstream>
using namespace std;

class Vehicle {
    public:
        string brand = "Ford";
        void honk() {
            cout << "Tuut, tuut\n";
        }
};

class Car: public Vehicle {
    public:
        string model = "Mustang";
};


class myClass {
    public:
        void myFunction() {
            cout << "Some content in parent class.\n";
        }
};

class myOtherClass {
    public:
        void myOtherFunction() {
            cout << "Some content in another class.\n";
        }
};

class myChild: public myClass, public myOtherClass {

};



class myGrandChild: public myChild {

};


class Employee {
    protected:
        int salary;
};

class Programmer: public Employee {
    public:
        int bonus;
        void setSalary(int s) {
            salary = s;
        }
        int getSalary() {
            return salary;
        }
};


class Animal {
    public:
        void animalSound() {
            cout << "the animal makes a sound\n";
        }
};

class Pig: public Animal {
    public:
        void animalSound() {
            cout << "the pig says: wee wee\n";
        }
};

class Dog: public Animal {
    public:
        void animalSound() {
            cout << "the dog says: bow bow\n";
        }
};

class BaseParams {
    public:
        double getTime() const {return _t;};
        void setTime(double t0) {_t = t0;};

        double _t = 0;
};
int main() {
    /* Car mycar; */
    /* mycar.honk(); */
    /* cout << mycar.brand + " " + mycar.model; */
    /* myChild myobj; */
    /* myobj.myFunction(); */
    /* myobj.myOtherFunction(); */
    /* Programmer guisong; */
    /* guisong.setSalary(1000); */
    /* guisong.bonus = 2000; */
    /* cout << "salary:" << guisong.getSalary() << endl; */
    /* cout << "salary:" << guisong.salary << endl; */
    /* cout << "bonus:" << guisong.bonus << endl; */
    /* Animal myAnimal; */
    /* Dog myDog; */
    /* Pig myPig; */
    
    /* myAnimal.animalSound(); */
    /* myDog.animalSound(); */
    /* myPig.animalSound(); */
    
    /* BaseParams myparam; */
    /* myparam._t = 10; */
    /* cout << myparam._t; */ 
    ofstream MyFile("filename.txt");
    MyFile << "Files can be tricky, but it is fun enough!";
    MyFile.close();
    return 0;
}


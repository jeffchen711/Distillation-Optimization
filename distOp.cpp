#include "distOp.h"
using namespace std;

//This program finds the optimum distillation values
//using the McCabe-Thiele method.
int main() {
	cout << "Welcome to the McCabe-Thiele Calculator!" << endl;
	cout << "Please enter alpha, zf, xb, and xd separate by lines." << endl;

	double a=0, zf=0, xb=0, xd=0;
	/*cin >> a;
	cin >> zf;
	cin >> xb;
	cin >> xd;
*/

    a=1.7, zf=0.6, xb=0.1, xd=0.8143;
    
    NodeArray *array = new NodeArray(250000, a, zf, xb, xd);
    
    array->calc();
    array->print();
    
    free(array);
    
    cout << "Process complete." << endl;
    return 0;
}
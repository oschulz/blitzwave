#include <iostream>
#include <time.h>


int main() {
	using namespace std;

	for (int e=7; e<=11; ++e) {
		for (int stride=(1<<e)-3; stride<=(1<<e)+3; ++stride) {
			int n = 345;
	
			int *data = new int[n*stride];

			double timer = clock();

			for (int dummy=0; dummy<100000; ++dummy) {
				for (int i=0; i<n; ++i) {
					data[stride*i] = 0;
				}
			}
			
			cout << "n: " << n << "\tstride: " << stride
			     << "\ttime: " << (clock() - timer)/100000 << endl;

			delete [] data;			
		}
		
		cout << endl;
	}
	
	return 0;
}

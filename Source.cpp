#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
using namespace std;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::nanoseconds;
using std::chrono::system_clock;

struct Arr {
	int* A;
	int length=0;
};

void HeapSort(Arr* A);
void InsertionSort(Arr* A);

void FillArrayRand(fstream *f) {
	for (long i = 0; i < 1000; i++) {
		std::time(nullptr);
		int r = rand();
		*f << r << "\n";
	}
}
void FillArrSort(fstream* f) {
	for (long i = 0; i < 1000; i++) {
		*f << i << "\n";
	}
}
void FillArrUnSort(fstream* f) {
	for (long i = 9000; i > 0; i--) {
		*f << i << "\n";
	}
}
void FillArrAlSort(fstream* f) {
	int n = 9000;
	for (long i = 0; i < n; i++) {
		if(i==n/2)*f << i+1 << "\n";
		else if(i==n/2 +1)*f << i -1<< "\n";
		else *f << i  << "\n";
	}
	
}

void IntroHeapify(int arr[], int n, int i) {
	int largest = i;
	int l = 2 * i + 1; 
	int r = 2 * i + 2;
	if (l < n && arr[l] >arr[largest])
		largest = l;
	if (r < n && arr[r] >arr[largest])
		largest = r;
	if (largest != i)
	{
		int tmp;
		tmp = arr[i];
		arr[i] = arr[largest];
		arr[largest] = tmp;
		IntroHeapify(arr, n, largest);
	}
}
void IntroHeapSort(int arr[], int size) {
	for (int i = size / 2 - 1; i >= 0; i--)
		IntroHeapify(arr, size, i);
	for (int i = size - 1; i >= 0; i--)
	{
		int tmp;
		tmp = arr[0];
		arr[0] = arr[i];
		arr[i] = tmp;
		IntroHeapify(arr, i, 0);
	}
}
void Fill(Arr *A, fstream* f) {
	int tmp; long length = 0;
	(*f).open("Array.txt", ios::in);
	(*f) >> tmp;
	while (!(*f).eof()) {
		length++; (*f) >> tmp;
	}
	(*A).A = new int[length];
	(*A).length = length;
	(*f).close();
	(*f).open("Array.txt", ios::in);
	length = 0;
	(*f) >> tmp;
	while (!(*f).eof()) {
		(*A).A[length] = tmp;
		(*f) >> tmp;
		length++;

	}
	(*f).close();
}
void Delete(Arr* A) {
	delete[](*A).A;
}
void OutputFile(Arr* A) {
	fstream f;
	f.open("Output Array.txt", ios::out);
	for (int i = 0; i < (*A).length; i++) {
		f << (*A).A[i] << endl;
	}
	f.close();
	
}
void Output(Arr* A) {
	cout << "\nCurrent array:\n";
	for (int i = 0; i < (*A).length; i++) {
		cout << (*A).A[i]<<" ";
	}
	cout << endl;
}
void Copy(Arr* A, Arr* B) {
	(*B).length = (*A).length;
	(*B).A = new int[(*B).length];
	for (int i = 0; i < (*B).length; i++)(*B).A[i] = (*A).A[i];

}
void Swap(Arr* A, int num1, int num2) {
	int tmp;
	tmp = (*A).A[num1];
	(*A).A[num1] = (*A).A[num2];
	(*A).A[num2] = tmp;
}
void Merge(Arr* A, Arr* B, Arr* C) {
	int i = 0; int j = 0; int n = 0;
	while (i < (*A).length  && j < (*B).length) {
		if ((*A).A[i] > (*B).A[j]) {
			(*C).A[n] = (*B).A[j];
			n++; j++;
		}
		else {
			(*C).A[n] = (*A).A[i];
			n++; i++;
		}
	}
	if (i != (*A).length) {
		while (i != (*A).length) {
			(*C).A[n] = (*A).A[i];
			n++; i++;
		}
	}
	else {
		while (j != (*B).length) {
			(*C).A[n] = (*B).A[j];
			n++; j++;
		}
	}
}
void Heapify(Arr* A, int n, int i) {
	int largest = i;
	int l = 2 * i + 1; 
	int r = 2 * i + 2; 
	if (l < n && (*A).A[l] >(*A).A[largest])
		largest = l;
	if (r < n && (*A).A[r] >(*A).A[largest])
		largest = r;
	if (largest != i)
	{
		swap((*A).A[i], (*A).A[largest]);
		Heapify(A, n, largest);
	}
}
void TimInsertionSort(Arr *A, int left, int right) {
	for (int i = left + 1; i <= right; i++) {
		int temp = (*A).A[i];
		int j = i - 1;
		while (j >= left && (*A).A[j] > temp) {
			(*A).A[j + 1] = (*A).A[j];
			j--;
		}
		(*A).A[j + 1] = temp;
	}
}
void TimMerge(Arr *A, int l, int m , int r) {
	int len1 = m - l + 1, len2 = r - m;
	int* left = new int[len1];int* right = new int[len2];
	for (int i = 0; i < len1; i++)
		left[i] = (*A).A[l + i];
	for (int i = 0; i < len2; i++)
		right[i] = (*A).A[m + 1 + i];

	int i = 0;
	int j = 0;
	int k = l;
	while (i < len1 && j < len2) {
		if (left[i] <= right[j]) {
			(*A).A[k] = left[i];
			i++;
		}
		else {
			(*A).A[k] = right[j];
			j++;
		}
		k++;
	}
	while (i < len1) {
		(*A).A[k] = left[i];
		k++;
		i++;
	}
	while (j < len2) {
		(*A).A[k] = right[j];
		k++;
		j++;
	}
}
int* Partition(int arr[], int low, int high)
{
	int pivot = arr[high];  
	int i = (low - 1);  
	for (int j = low; j <= high - 1; j++)
	{
		if (arr[j] <= pivot)
		{
			i++;
			swap(arr[i], arr[j]);
		}
	}
	swap(arr[i + 1], arr[high]);
	return (arr + i + 1);
}
int* MedianOfThree(int* a, int* b, int* c)
{
	
		
			if (*a < *b && *b < *c)
				return (b);

			if (*a < *c && *c <= *b)
				return (c);

			if (*b <= *a && *a < *c)
				return (a);

			if (*b < *c && *c <= *a)
				return (c);

			if (*c <= *a && *a < *b)
				return (a);

			if (*c <= *b && *b <= *a)
				return (b);
		
	
	
}
void IntroUtil(Arr *A, int* begin, int* end, int depthLimit) {
	int size = end - begin;
	if (size < 16)
	{
		InsertionSort(A);
		return;
	}
	if (depthLimit == 0)
	{
		Heapify(A, *begin, *end + 1);
		HeapSort(A);
		return;
	}
	if (*begin >= 0 && *end>*begin) {
		int* pivot = MedianOfThree(begin, begin + size / 2, end);
		swap(pivot, end);
	

	int pivot1 =  (*A).A[*end];    
	int i = *(begin) + size/2 - 1;  

	for (int j = *(begin)+size / 2; j <= *end - 1; j++)
	{
		if ((*A).A[j] <= pivot1)
		{
			i++;
			swap((*A).A[i], (*A).A[j]);
		}
	}
	swap((*A).A[i + 1], (*A).A[*end]);
	int* partitionPoint = (*A).A + i + 1;
	IntroUtil(A, begin, partitionPoint-1 , depthLimit - 1);
	IntroUtil(A, partitionPoint+1 , end, depthLimit - 1);
	}
}

void SelectionSort(Arr* A) {

	int min; int num; bool flag;
	for (int i =0; i < (*A).length; i++) {
		min = (*A).A[i];
		flag =0;
		for (int j = i+1; j < (*A).length; j++) {
			if ((*A).A[j] < min) {
				min = (*A).A[j];
				num = j;
				flag = 1;
			}
		}
		if(flag==1)Swap(A, i, num);
	}
}
void InsertionSort(Arr* A) {
	for (int i = 0; i < (*A).length; i++) {
		for (int j = i; j>0 && (*A).A[j - 1]>(*A).A[j]; j--) {
			Swap(A, j - 1, j);
		}
	}
}
void BubbleSort(Arr* A) {
	for (int i = 0; i < (*A).length-1; i++) {
		for (int j = 0; j < (*A).length-1; j++) {
			if ((*A).A[j] > (*A).A[j + 1]) {
				Swap(A, j, j + 1);
			}
		}
	}
}
void MergeSort(Arr* A,int leng) {
	if (leng == 1) {  }
	else {
		Arr left, right;
		left.length = leng / 2;
		left.A = new int[left.length];
		right.length = leng / 2;
		if (leng % 2 != 0)right.length++;
		right.A = new int[right.length];
		for (int i = 0; i < left.length; i++)left.A[i] = (*A).A[i];
		for (int i = left.length; i < (*A).length; i++)right.A[i- left.length] = (*A).A[ i];
		MergeSort(&left, left.length);
		MergeSort(&right, right.length);
		Merge(&left, &right, A);
	}
}
void QuickSort(Arr* A, int start, int end) {
	if (start >= end)
		return;

	int pivot = (*A).A[start];

	int count = 0;
	for (int i = start + 1; i <= end; i++) {
		if ((*A).A[i] <= pivot)
			count++;
	}

	int pivotIndex = start + count;
	swap((*A).A[pivotIndex], (*A).A[start]);

	int i = start, j = end;

	while (i < pivotIndex && j > pivotIndex) {

		while ((*A).A[i] <= pivot) {
			i++;
		}

		while ((*A).A[j] > pivot) {
			j--;
		}

		if (i < pivotIndex && j > pivotIndex) {
			swap((*A).A[i++], (*A).A[j--]);
		}
	}

	int p = pivotIndex;

	QuickSort(A, start, p - 1);

	QuickSort(A, p + 1, end);
}
void ShellSort(Arr* A) {
	int i, j, step, m = 0, n=0;
	int tmp;
	//for (m = (log2((*A).length)); m > 0; m--) {
		//	for (step = (*A).length / 2; step > 0; step /= 2)
		//step = pow(2, m) - 1;
	for (step = (*A).length / 2; step > 0; step /= 2)
			for (i = step; i < (*A).length; i++)
			{
				tmp = (*A).A[i];
				for (j = i; j >= step; j -= step)
				{
					if (tmp < (*A).A[j - step])
						(*A).A[j] = (*A).A[j - step];
					else
						break;
				}
				(*A).A[j] = tmp;
			}
		
	
	//}
}
void HeapSort(Arr* A) {
	for (int i = (*A).length / 2 - 1; i >= 0; i--)
		Heapify(A,(*A).length, i);
	for (int i = (*A).length - 1; i >= 0; i--)
	{
		Swap(A, 0, i);
		Heapify(A, i, 0);
	}
}
void TimSort(Arr* A) {
	int RUN = 32; int min;
	for (int i = 0; i < (*A).length; i += RUN) {
		if (i + RUN - 1 < (*A).length - 1) min = i + RUN - 1;
		else min = (*A).length - 1;
		TimInsertionSort(A, i, min);
	}
	for (int size = RUN; size < (*A).length; size = 2 * size) {
		for (int left = 0; left < (*A).length; left += 2 * size) {
			int mid = left + size - 1;
			int right;
			if (left + 2 * size - 1 < (*A).length - 1) right = left + 2 * size - 1;
			else right = (*A).length - 1;
			if (mid < right)
				TimMerge(A, left, mid, right);
		}
	}
}
void IntroSort(Arr* A, int* begin, int* end) {
	int depthLimit = 2 * log(end - begin) - 1;
	IntroUtil(A, begin, end, depthLimit);

}

void Menu(Arr* A, Arr* B) {
	char c; bool exit=0;
	fstream h;
	h.open("Output time.txt", ios::app);
	for (;;) {
		Copy(A, B);
		cout << "Please choose sort:\n1) Selection sort\n2) Insertion sort\n3) Bubble sort \n4) Merge sort\n5) Quick sort\n6) Shell sort\n7) Heap sort\n8) Tim sort\n9) Intro sort\n0) Exit.\n" ;
		cin >> c;
		uint64_t ms = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
		switch (c) {
		case '1':
			SelectionSort(B);
			break;
		
		case '2':
			InsertionSort(B);
			break;
		case '3':
			BubbleSort(B);
			break;
		case '4':
			MergeSort(B, (*B).length);
			break;
		case '5':
			QuickSort(B, 0, (*B).length-1);
			break;
		case '6':
			ShellSort(B);
			break;
		case '7':
			HeapSort(B);
			break;
		case '8':
			TimSort(B);
			break;
		case '9':
			IntroSort(B,0, &B->length-1);
			break;
		case '0':
			cout << "\nBye!\n";
			exit = 1;
			break;
		default: 
			cout << "\nWrong input.\n";
			break;
		}
		if (exit == 1)break;
		uint64_t m = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
		
		
		//cout << "\nElapsed time: " << m - ms << " ns.\n";
		h <<  m - ms << " ";
		OutputFile(B);
		//Output(A);
		//Output(B);
		
		
	}
	h << endl;
	h.close();
}

int main() {
	fstream f,g,p;
	Arr A, B;
	g.open("Output Array.txt", ios::out);
	f.open("Array.txt", ios::out);
	p.open("Output time.txt", ios::out);
	if (!f.is_open() || !g.is_open())cout << "Could not open files. Please check for \"file Array.txt\" and \"Output Array.txt\"\n";
	else {
		FillArrayRand(&f);
		//FillArrSort(&f);
		//FillArrAlSort(&f);

		//FillArrUnSort(&f);

			f.close();
			g.close();
			cout << "Hello! Filling Array...\n";
			Fill(&A, &f);
			//Output(&A);
			Menu(&A, &B);
			Delete(&A);
			Delete(&B);
	}
}
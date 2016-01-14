#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

//#include <stdio.h>
//#include <stdlib.h>
#include <time.h>

struct Header{
	char chunkID[4];
	int chunkSize;
	char format[4];

	char subChunk1ID[4];
	int subChunk1Size;
	short audioFormat;
	short numChannels;
	int sampleRate;
	int byteRate;
	short blockAlign;
	short bitsPerSample;

	int channelSize;
	char subChunk2ID[4];
	int subChunk2Size;
};

struct Header dryHeader;
struct Header irHeader;

double* data;
double* irdata;
double* outdata;

int dryNumSamples;
int irNumSamples;
int outNumSamples;



//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;

	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j + 1];
				tempi = wr * data[j + 1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}


void convolve()
{
	 double * x;
	 double * h;
	 double * y;
	 
	int nn = (int)pow(2, (int) log2(dryNumSamples) + 1);
	int nn_2 = 2 * nn; 
	
	 
	 h = (double*) malloc(sizeof(double) * nn_2); 
	 x = (double*) malloc(sizeof(double) * nn_2); 
	 y = (double*) malloc(sizeof(double) * nn_2); 
	 outdata = (double*)malloc(sizeof(double) * nn);

	 
	 // initialize and zero pad the arrays
	for(int i = 0; i < nn_2; i++)
	{
		x[i] = 0;
	}

	for(int i = 0; i < nn_2; i++)
	{
		h[i] = 0;
	}

	for(int i = 0; i < nn_2; i++)
	{
		y[i] = 0;
	}

	
	//write the dry data.
	for(int i = 0; i < dryNumSamples; i++)
	{
		x[2 * i] = data[i];
	}

	//write the ir data
	for(int i = 0; i < irNumSamples; i++)
	{
		h[2 * i] = irdata[i];
	}

	
	four1(h - 1, nn, 1);
	four1(x - 1, nn, 1);
	 
	double MAX_VAL = 32767;
	// Complex multiplication i think (?)
	for(int i = 0 ; i < nn_2; i+=2)
	{
		y[i] = (x[i]/MAX_VAL * h[i]/MAX_VAL) - (x[i + 1]/MAX_VAL * h[i + 1]/MAX_VAL);
		y[i + 1] = (x[i + 1]/MAX_VAL * h[i]/MAX_VAL) + (x[i]/MAX_VAL * h[i + 1]/MAX_VAL);
	}


	// inverse
	four1(y - 1, nn, -1);

	//have to put the shit back in the right spots
	for(int i = 0; i < nn_2; i++)
	{
		outdata[i] = y[i*2]/(nn_2 * 2);
	}

	
	free(x);
	free(h);
}



/////////////////////////////////////////////////////////////////////////////////////////
void print()
{
	printf("\n=============IN HEADER INFO =============\n");
	printf(" chunkID:%s\n", dryHeader.chunkID);
	printf(" chunkSize:%d\n", dryHeader.chunkSize);
	printf(" format:%s\n", dryHeader.format);
	printf(" subChunk1ID:%s\n", dryHeader.subChunk1ID);
	printf(" subChunk1Size:%d\n", dryHeader.subChunk1Size);
	printf(" audioFormat:%d\n", dryHeader.audioFormat);
	printf(" numChannels:%d\n", dryHeader.numChannels);
	printf(" sampleRate:%d\n", dryHeader.sampleRate);
	printf(" byteRate:%d\n", dryHeader.byteRate);
	printf(" blockAlign:%d\n", dryHeader.blockAlign);
	printf(" bitsPerSample:%d\n", dryHeader.bitsPerSample);
	printf(" subChunk2ID:%s\n", dryHeader.subChunk2ID);
	printf(" subChunk2Size:%d\n", dryHeader.subChunk2Size);
}

void printIR()
{
	printf("\n=============IR HEADER INFO =============\n");
	printf(" chunkID:%s\n", irHeader.chunkID);
	printf(" chunkSize:%d\n", irHeader.chunkSize);
	printf(" format:%s\n", irHeader.format);
	printf(" subChunk1ID:%s\n", irHeader.subChunk1ID);
	printf(" subChunk1Size:%d\n", irHeader.subChunk1Size);
	printf(" audioFormat:%d\n", irHeader.audioFormat);
	printf(" numChannels:%d\n", irHeader.numChannels);
	printf(" sampleRate:%d\n", irHeader.sampleRate);
	printf(" byteRate:%d\n", irHeader.byteRate);
	printf(" blockAlign:%d\n", irHeader.blockAlign);
	printf(" bitsPerSample:%d\n", irHeader.bitsPerSample);
	printf(" subChunk2ID:%s\n", irHeader.subChunk2ID);
	printf(" subChunk2Size:%d\n", irHeader.subChunk2Size);
}


int loadWave(char* filename)
{
	FILE* in = fopen(filename, "rb");

	if (in != NULL)
	{		
		printf("Reading %s...\n",filename);

		fread(dryHeader.chunkID, 1, 4, in);
		fread(&dryHeader.chunkSize, 1, 4, in);
		fread(dryHeader.format, 1, 4, in);

		//sub chunk 1
		fread(dryHeader.subChunk1ID, 1, 4, in);
		fread(&dryHeader.subChunk1Size, 1, 4, in);
		fread(&dryHeader.audioFormat, 1, 2, in);
		fread(&dryHeader.numChannels, 1, 2, in);
		fread(&dryHeader.sampleRate, 1, 4, in);
		fread(&dryHeader.byteRate, 1, 4, in);
		fread(&dryHeader.blockAlign, 1, 2, in);
		fread(&dryHeader.bitsPerSample, 1, 2, in);		
		
		//read extra bytes
		if(dryHeader.subChunk1Size == 18)
		{
			short empty;
			fread(&empty, 1, 2, in);		
		}
		
		//sub chunk 2
		fread(dryHeader.subChunk2ID, 1, 4, in);
		fread(&dryHeader.subChunk2Size, 1, 4, in);

		//read data		
		int bytesPerSample = dryHeader.bitsPerSample/8;
		dryNumSamples = dryHeader.subChunk2Size / bytesPerSample;
		data = (double*) malloc(sizeof(double) * dryNumSamples);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		double MAX_VAL = 32767;
		int i=0;
		short sample=0;
		while(fread(&sample, 1, bytesPerSample, in) == bytesPerSample)
		{		
			data[i++] = (double)sample;
			sample = 0;			
		}
		
		fclose(in);
		printf("Closing %s...\n",filename);
	}
	else
	{
		printf("Can't open file\n");
		return 0;
	}
	return 1;
}

int loadIRWave(char* filename)
{
	FILE* in = fopen(filename, "rb");

	if (in != NULL)
	{		
		printf("Reading %s...\n",filename);

		fread(irHeader.chunkID, 1, 4, in);
		fread(&irHeader.chunkSize, 1, 4, in);
		fread(irHeader.format, 1, 4, in);

		//sub chunk 1
		fread(irHeader.subChunk1ID, 1, 4, in);
		fread(&irHeader.subChunk1Size, 1, 4, in);
		fread(&irHeader.audioFormat, 1, 2, in);
		fread(&irHeader.numChannels, 1, 2, in);
		fread(&irHeader.sampleRate, 1, 4, in);
		fread(&irHeader.byteRate, 1, 4, in);
		fread(&irHeader.blockAlign, 1, 2, in);
		fread(&irHeader.bitsPerSample, 1, 2, in);		
		
		//read extra bytes
		if(irHeader.subChunk1Size == 18)
		{
			short empty;
			fread(&empty, 1, 2, in);		
		}
		
		//sub chunk 2
		fread(irHeader.subChunk2ID, 1, 4, in);
		fread(&irHeader.subChunk2Size, 1, 4, in);

		//read data		
		int bytesPerSample = irHeader.bitsPerSample/8;
		irNumSamples = irHeader.subChunk2Size / bytesPerSample;
		irdata = (double*) malloc(sizeof(double) * irNumSamples);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		double MAX_VAL = 32767;
		int i=0;
		short sample=0;
		while(fread(&sample, 1, bytesPerSample, in) == bytesPerSample)
		{		
			irdata[i++] = (double)sample;
			sample = 0;			
		}
		
		fclose(in);
		printf("Closing %s...\n",filename);
	}
	else
	{
		printf("Can't open file\n");
		return 0;
	}
	return 1;
}

int saveWave(char* filename)
{
	FILE* out = fopen(filename, "wb");

	if (out != NULL)
	{		
		printf("\nWriting %s...\n",filename);

		fwrite(dryHeader.chunkID, 1, 4, out);
		fwrite(&dryHeader.chunkSize, 1, 4, out);
		fwrite(dryHeader.format, 1, 4, out);

		//sub chunk 1
		fwrite(dryHeader.subChunk1ID, 1, 4, out);
		fwrite(&dryHeader.subChunk1Size, 1, 4, out);
		fwrite(&dryHeader.audioFormat, 1, 2, out);
		fwrite(&dryHeader.numChannels, 1, 2, out);
		fwrite(&dryHeader.sampleRate, 1, 4, out);
		fwrite(&dryHeader.byteRate, 1, 4, out);
		fwrite(&dryHeader.blockAlign, 1, 2, out);
		fwrite(&dryHeader.bitsPerSample, 1, 2, out);		
		
		//read extra bytes
		if(dryHeader.subChunk1Size == 18)
		{
			short empty = 0;
			fwrite(&empty, 1, 2, out);		
		}
		
		//sub chunk 2
		fwrite(dryHeader.subChunk2ID, 1, 4, out);
		fwrite(&dryHeader.subChunk2Size, 1, 4, out);

		//read data		
		int bytesPerSample = dryHeader.bitsPerSample / 8;
		int sampleCount =  dryHeader.subChunk2Size / bytesPerSample;

		double maxSample = -1;
		double MAX_VAL = 32767;	//FIXME: find based on bits per sample
			
		for(int i=0; i < dryNumSamples; i++)
		{			

			 if(i==0)
			 	maxSample = outdata[0];
			 else if(outdata[i] > maxSample)
			 	maxSample = outdata[i];
		}		
		
		//scale and re write the data
		for(int i=0; i<dryNumSamples; i++)
		{
			outdata[i] = (outdata[i] / maxSample) ;
			short sample = (short) (outdata[i] * MAX_VAL);
			fwrite(&sample, 1, bytesPerSample, out);
		}
		
		//clean up
		//	free(newData);
		fclose(out);
		printf("Closing %s...\n",filename);
	}
	else
	{
		printf("Can't open file\n");
		return 0;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	printf("=============================================================\n");
	printf("=============================================================\n");
	char *inputFileName;
	char *irFileName;
	char *outputFileName;

	if(argc == 4)
	{
		printf("Correct number of inputs\n");
	}
	else
	{
		printf("Usage: ./convolve <Dry recording> <Impulse Response> <Output filename\n");
		return 0;
	}

	inputFileName = argv[1];
	irFileName = argv[2];
	outputFileName = argv[3];

	printf("Input file is: %s\n", inputFileName);
	printf("Impulse response file is: %s\n", irFileName);
	printf("Output file is: %s\n", outputFileName);

	clock_t start = clock();
	
	if(loadWave(inputFileName))
		print();
		
	if(loadIRWave(irFileName));
		printIR();

	outNumSamples = dryNumSamples + irNumSamples - 1;
	convolve();
		
	saveWave(outputFileName);

	clock_t end = clock();


	free(data);
	free(irdata);
	free(outdata);

	

	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("Time to run entire program in seconds: %f\n", seconds);
}
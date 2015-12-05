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

short* data;
short* irdata;

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



// Creates a sine tone with the specified harmonic number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0).

void createComplexSine(double data[], int size, int harmonicNumber)
{
	int i, ii;

	for (i = 0, ii = 0; i < size; i++, ii += 2) {
		data[ii] = sin((double)harmonicNumber * (double)i * TWO_PI / (double)size);
		data[ii + 1] = 0.0;
	}
}



// Creates a cosine tone with the specified harmonic number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0).

void createComplexCosine(double data[], int size, int harmonicNumber)
{
	int i, ii;

	for (i = 0, ii = 0; i < size; i++, ii += 2) {
		data[ii] = cos((double)harmonicNumber * (double)i * TWO_PI / (double)size);
		data[ii + 1] = 0.0;
	}
}



// Creates a sawtooth wave, where each harmonic has
// the amplitude of 1 / harmonic_number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0)

void createComplexSawtooth(double data[], int size)
{
	int i, ii, j;

	//  Calculate waveform using additive synthesis
	for (i = 0, ii = 0; i < size; i++, ii += 2) {
		data[ii] = 0.0;
		data[ii + 1] = 0.0;
		for (j = 1; j <= size / 2; j++) {
			data[ii] +=
				(cos((double)j * (double)i * TWO_PI / (double)size)) / (double)j;
		}
	}
}



// Display the real and imaginary parts
// the data contained in the array.

void displayComplex(double data[], int size)
{
	int i, ii;

	printf("\t\tReal part \tImaginary Part\n");

	for (i = 0, ii = 0; i < size; i++, ii += 2)
		printf("data[%-d]: \t%.6f \t%.6f\n", i, data[ii], data[ii + 1]);

	printf("\n");
}



// Performs the DFT on the input data,
// which is assumed to be a real signal.
// That is, only data at even indices is
// used to calculate the spectrum.

void complexDFT(double x[], int N)
{
	int n, k, nn;
	double omega = TWO_PI / (double)N;
	double *a, *b;

	// Allocate temporary arrays
	a = (double *)calloc(N, sizeof(double));
	b = (double *)calloc(N, sizeof(double));

	// Perform the DFT
	for (k = 0; k < N; k++) {
		a[k] = b[k] = 0.0;
		for (n = 0, nn = 0; n < N; n++, nn += 2) {
			a[k] += (x[nn] * cos(omega * n * k));
			b[k] -= (x[nn] * sin(omega * n * k));
		}
	}

	// Pack result back into input data array
	for (n = 0, k = 0; n < N * 2; n += 2, k++) {
		x[n] = a[k];
		x[n + 1] = b[k];
	}

	// Free up memory used for arrays
	free(a);
	free(b);
}



// Takes the results from a DFT or FFT, and
// calculates and displays the amplitudes of
// the harmonics.

void postProcessComplex(double x[], int N)
{
	int i, k, j;
	double *amplitude, *result;

	// Allocate temporary arrays
	amplitude = (double *)calloc(N, sizeof(double));
	result = (double *)calloc(N, sizeof(double));

	// Calculate amplitude
	for (k = 0, i = 0; k < N; k++, i += 2) {
		// Scale results by N
		double real = x[i] / (double)N;
		double imag = x[i + 1] / (double)N;
		// Calculate amplitude
		amplitude[k] = sqrt(real * real + imag * imag);
	}

	// Combine amplitudes of positive and negative frequencies
	result[0] = amplitude[0];
	result[N / 2] = amplitude[N / 2];
	for (k = 1, j = N - 1; k < N / 2; k++, j--)
		result[k] = amplitude[k] + amplitude[j];


	// Print out final result
	printf("Harmonic \tAmplitude\n");
	printf("DC \t\t%.6f\n", result[0]);
	for (k = 1; k <= N / 2; k++)
		printf("%-d \t\t%.6f\n", k, result[k]);
	printf("\n");

	// Free up memory used for arrays
	free(amplitude);
	free(result);
}



int main()
{
	int i;
	double complexData[SIZE * 2];

	// Try the DFT on a sawtooth waveform
	createComplexSawtooth(complexData, SIZE);
	displayComplex(complexData, SIZE);
	complexDFT(complexData, SIZE);
	postProcessComplex(complexData, SIZE);

	// Try the FFT on the same data
	createComplexSawtooth(complexData, SIZE);
	displayComplex(complexData, SIZE);
	four1(complexData - 1, SIZE, 1);
	postProcessComplex(complexData, SIZE);
}


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
		data = (short*) malloc(sizeof(short) * dryNumSamples);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		
		int i=0;
		short sample=0;
		while(fread(&sample, 1, bytesPerSample, in) == bytesPerSample)
		{		
			data[i++] = sample;
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
		irdata = (short*) malloc(sizeof(short) * irNumSamples);
		
		//fread(data, 1, bytesPerSample*numSamples, in);
		
		int i=0;
		short sample=0;
		while(fread(&sample, 1, bytesPerSample, in) == bytesPerSample)
		{		
			irdata[i++] = sample;
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
		
		//impulse response - echo
		// int IRSize = 6;
		// float IR[IRSize];
		// IR[0] = 1.0;
		// IR[1] = 1.0;
		// IR[2] = 1.0;
		// IR[3] = 1.0;
		// IR[4] = 1.0;
		// IR[5] = 1.0;
		
		//write the data
		float* newData = (float*) malloc(sizeof(float) * (outNumSamples));// + IRSize - 1));
		float maxSample = -1;
		float MAX_VAL = 32767.f;	//FIXME: find based on bits per sample
			
		for(int i=0; i < dryNumSamples; ++i)
		{			

			//convolve
			for(int j=0; j < irNumSamples; ++j)
				newData[i+j] += ((float)data[i] / MAX_VAL) * ((float)irdata[j] / MAX_VAL);
			
			//Keep track of max value for scaling
			 if(i==0)
			 	maxSample = newData[0];
			 else if(newData[i] > maxSample)
			 	maxSample = newData[i];
		}		
		
		//scale and re write the data
		for(int i=0; i<outNumSamples; ++i)
		{
			newData[i] = (newData[i] / maxSample) ;
			short sample = (short) (newData[i] * MAX_VAL);
			fwrite(&sample, 1, bytesPerSample, out);
		}
		
		//clean up
		free(newData);
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

	clock_t start = clock();

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
	
	if(loadWave(inputFileName))
		print();
		
	if(loadIRWave(irFileName));
		printIR();

	outNumSamples = dryNumSamples + irNumSamples - 1;
		
	saveWave(outputFileName);
	free(data);
	free(irdata);

	clock_t end = clock();

	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("%f\n", seconds);
}
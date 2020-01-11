#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#define MaxRows 500
#define MaxCols 500

#define MaxMaskSize 10

using namespace std;

struct Mask {
public:
	double MaskArray[MaxMaskSize][MaxMaskSize];
	int Rows, Cols;
};

struct Image {
private:
	unsigned short SourceImage[MaxRows][MaxCols];
	int Rows, Cols;
	bool Available;

	void init() {
		Rows = 0;
		Cols = 0;
		Available = false;
	}
public:
	Image() {
		init();
	}

	Image(string Filename) {
		loadImage(Filename);
	}

	int SetPixel(unsigned short Value, int R, int C) {
		if (R > -1 && C > -1 && R < MaxRows && C < MaxCols) {
			SourceImage[R][C] = Value;
			return Value;
		}
		return -1;
	}

	bool isAvailable() {
		return Available;
	}

	int loadImage(string FileName) {
		init();
		// Code to load an image from a storage into array named Image in main memory
		ifstream input;		//handle to get info from inputfilestream
		input.open(FileName.c_str());
		if (!input)
			return -1;
		string format, comments; int maxcolor;
		do
		{
			getline(input, format, '\n');
			cout << format << endl;

			getline(input, comments, '\n');
			cout << comments << endl;
			input >> Cols;
			input >> Rows;
			cout << Rows << " " << Cols << endl;
			input >> maxcolor;
			cout << maxcolor;
			for (int i = 0; i < Rows; i++)
			{
				for (int j = 0; j <Cols; j++)
				{
					input >> SourceImage[i][j];
				}
			}

		}while (!input.eof());
		Available = true;
		return 0;
	}

	int saveImage(string FileName) {
		if (!Available)
			return 1;
		// Code to store the SourceImage or ProcessedImage into a file on some backing storage
		//  WhichImage 0 means copy source image and otherwise ProcessedImage
		ofstream writeout;
		writeout.open(FileName.c_str());
		writeout << "P2" << endl;
		writeout << "# created by Mubasharians" << endl;
		writeout << Cols << " " << Rows << endl;
		writeout << 255 << endl;
		for (size_t i = 0; i < Rows; i++)
		{
			for (size_t j = 0; j <Cols; j++)
			{
                writeout << SourceImage[i][j] << " ";
			}
			writeout << endl;
		}
		return 0;
	}

	void RotateImageAboutPoint(Image &Result, double ByAngle = 90, int Cx = 0, int Cy = 0) {
		// Code to compute a rotated version of original image
		// The image will be rotated about the center (Cx, Cy)
		// and the resultant image will be stored in Result
		ByAngle = (ByAngle / 180)*3.142;

		int xDash;
		int yDash;

		for (int i = 0; i < Rows; i++)
		{
			for (int j = 0; j < Cols; j++)
			{
				xDash = (j-Cx)*cos(ByAngle) - (i-Cy)*sin(ByAngle) + 0.5 + Cx;
				yDash = (j-Cx)*sin(ByAngle) + (i-Cy)*cos(ByAngle) + 0.5 + Cy;

				if (xDash > 0 && xDash <= Cols && yDash > 0 && yDash <= Rows)
				{
					Result.SourceImage[i][j] = SourceImage[xDash][yDash];
				}
			}
		}
        Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;

}

	void TranslateImage(Image &Result, int Tx, int Ty) {
		// Code to compute a translated version of original image
		// The image will be translated by Tx and Ty respectively along X and Y axes
		int xDash;
		int yDash;
        Cols = Cols + Ty ;
        Rows = Rows + Tx;
		for (int i = 0; i < Rows ; i++)
		{
			for (int j = 0; j < Cols ; j++)
			{
				xDash = (i+Tx) ;
				yDash = (j+Ty) ;
				Result.SourceImage[i][j] = SourceImage[xDash][yDash];
			}
		}
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.Available = true;


	}

	void ScaleImage(Image &Result, double SX, double SY) {
		// Code to compute a scaled version of the original image
		// The image will be scaled by ratios Sx and Sy respectively along X and Y axes
		for (double i = 0, x = 0; i < Rows;i += 1/SY, x++)
            for (double j = 0, y = 0; j < Cols;j += 1/SX, y++)
                Result.SourceImage[(int)x][(int)y] = SourceImage[(int)i][(int)j];

        Result.Cols = Cols * SX;
        Result.Rows = Rows * SY;
        Result.Available = true;

	}

	void LinearTransform(Image &Result, double TransformationMatrix[2][2]) {
		// Code to compute a transformed version of the original image
		// The transformation matrix will be specified in LinearTransformation matrix
		int xDash = 0;
		int yDash = 0;
		for(int i = 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                xDash = TransformationMatrix[0][0]*i + TransformationMatrix[0][1]*j;
                yDash = TransformationMatrix[1][0]*i + TransformationMatrix[1][1]*j;
                Result.SourceImage[xDash][yDash] = SourceImage[i][j];
            }
        }
        Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;

	}

	void FlipImage(Image& Result, int Axis) {
		// Code to compute a flipped version of the original image
		// The Axis parameter will specify the flip axis.
		// 1 means flip along X-axis and otherwise flip along Y-axis
        if(0 == Axis)
        {
            for(int r = 0 ; r < (Rows);r++)
            {
                for(int c = 0 ; c != Cols ; c++)
                {
                    Result.SourceImage[Rows-1-r][c] = SourceImage[r][c];
                }
            }
        }
        else
        {
            for(int r = 0 ; r < Rows;r++)
            {
                for(int c = 0 ; c < Cols ; c++)
                {
                    Result.SourceImage[r][Cols-c-1] = SourceImage[r][c];
                }
            }
        }
        Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;
	}

	void QuantizeImage(Image &Result , int RangeArray[][2]) {

        for(int i = 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                if( SourceImage[i][j] <= RangeArray[0][1] && SourceImage[i][j] > 0)
                    Result.SourceImage[i][j] = (RangeArray[0][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[1][1] && SourceImage[i][j] > RangeArray[0][1])
                    Result.SourceImage[i][j] = (RangeArray[0][1] + RangeArray[1][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[2][1] && SourceImage[i][j] > RangeArray[1][1])
                    Result.SourceImage[i][j] = (RangeArray[2][1] + RangeArray[1][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[3][1] && SourceImage[i][j] > RangeArray[2][1])
                    Result.SourceImage[i][j] = (RangeArray[3][1] + RangeArray[2][1] ) / 2;
                else if( SourceImage[i][j] <= RangeArray[4][1] && SourceImage[i][j] > RangeArray[3][1])
                    Result.SourceImage[i][j] = (RangeArray[4][1] + RangeArray[3][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[5][1] && SourceImage[i][j] > RangeArray[4][1])
                    Result.SourceImage[i][j] = (RangeArray[5][1] + RangeArray[4][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[6][1] && SourceImage[i][j] > RangeArray[5][1])
                    Result.SourceImage[i][j] = (RangeArray[6][1] + RangeArray[5][1]) / 2;
                else if( SourceImage[i][j] <= RangeArray[7][1] && SourceImage[i][j] > RangeArray[6][1])
                    Result.SourceImage[i][j] = (RangeArray[7][1] + RangeArray[6][1]) / 2;
            }
        }
        Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;
	}

	void MeanFilterImage(Image &Result, int Neighborhood = 3) {
        for(int i = 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                Result.SourceImage[i][j] = ( SourceImage[i][j] + SourceImage[i-1][j+1] + SourceImage[i-1][j]
                                            + SourceImage[i-1][j-1] + SourceImage[i][j+1] + SourceImage[i][j-1]
                                            + SourceImage[i+1][j+1] + SourceImage[i+1][j-1] + SourceImage[i+1][j])/9;
            }
        }
        Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;
	}

	void MedianFilterImage(Image &Result, int Neighborhood = 3) {
        int N[9];
        for(int i= 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                N[1] = SourceImage[i][j];
                N[2] = SourceImage[i-1][j-1];
                N[3] = SourceImage[i-1][j];
                N[4] = SourceImage[i-1][j+1];
                N[5] = SourceImage[i][j-1];
                N[6] = SourceImage[i][j+1];
                N[7] = SourceImage[i+1][j-1];
                N[8] = SourceImage[i+1][j];
                N[9] = SourceImage[i+1][j+1];
                for(int k = 1 ; k <=9 ; k++)
                {
                    for(int l = k+1 ; l <= 9 ; l++)
                    {
                        if(N[k] > N[l])
                        {
                            int temp = N[k];
                            N[k] = N[l];
                            N[l] = temp;
                        }
                    }
                }
                Result.SourceImage[i][j] = N[5];
            }
        }
	    Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;
	}

	void ComputeDerivativeImage(Image &Result, Mask DerivativeMask) {

        int sum = 0;
        for(int i= 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                sum += SourceImage[i][j] * 0;
                sum += SourceImage[i-1][j-1] * -1;
                sum += SourceImage[i-1][j] * -1;
                sum += SourceImage[i-1][j+1] * -1;
                sum += SourceImage[i][j-1] * -1;
                sum += SourceImage[i][j+1] * 1;
                sum += SourceImage[i+1][j-1] * 1;
                sum += SourceImage[i+1][j] * 1;
                sum += SourceImage[i+1][j+1] * 1;
                Result.SourceImage[i][j] = abs(sum);
                sum = 0;
            }
        }
	    Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;
	}

	void ApplyFilter(Image &Result, Mask FilterMask) {

	    int sum = 0;
        for(int i= 0 ; i < Rows ; i++)
        {
            for(int j = 0 ; j < Cols ; j++)
            {
                sum += SourceImage[i][j] * (FilterMask.MaskArray[1][1] * 100);
                sum += SourceImage[i-1][j-1] * (FilterMask.MaskArray[0][0] * 100);
                sum += SourceImage[i-1][j] * (FilterMask.MaskArray[0][1] * 100);
                sum += SourceImage[i-1][j+1] * (FilterMask.MaskArray[0][2] * 100);
                sum += SourceImage[i][j-1] * (FilterMask.MaskArray[1][0] * 100);
                sum += SourceImage[i][j+1] * (FilterMask.MaskArray[1][2] * 100);
                sum += SourceImage[i+1][j-1] * (FilterMask.MaskArray[2][0] * 100);
                sum += SourceImage[i+1][j] * (FilterMask.MaskArray[2][1] * 100);
                sum += SourceImage[i+1][j+1] * (FilterMask.MaskArray[2][2] * 100);
                Result.SourceImage[i][j] = abs(sum/100);
                sum = 0;
            }
        }

	    Result.Available = true;
	    Result.Cols = Cols;
	    Result.Rows = Rows;

	}
};

struct Menu {
private:
	string Options[30];
	int MaxOptions;
public:

	int GetMaxOptions() {
		return MaxOptions;
	}

	Menu(string FileName) {
		MaxOptions = 0;

		ifstream OptionsFile(FileName.c_str());

		if (!OptionsFile) {
			return;
		}
		while (!OptionsFile.eof()) {
			getline(OptionsFile, Options[MaxOptions++]);
		}
		return;
	}

	int ShowMenuAndGetChoice() {
		int Choice = 1;

		cout << "\nWel-Come To Mini-Project About Image Processing (Arrays)\n";
		do
		{
			if (Choice < 1 || Choice > MaxOptions) {
				cout << endl << " Please Make a Valid Choice\n\n Press enter to continue";
			}

			for (int i = 0; i< MaxOptions; i++)
				cout << endl << i + 1 << ":\t" << Options[i];

			cout << endl << endl << MaxOptions + 1 << ":\t Exit";

			cout << endl << "\nMake a choice from the Menu by Specifying Choice No ";
			cin >> Choice;
		} while (Choice < 1 || Choice > MaxOptions + 1);

		return Choice;
	}
};



void LOAD(Image &I) {
	string FileName = "";
	cout << "Specify Image File Name: ";
	do {
		getline(cin, FileName);
	} while (FileName.length() < 1);

	if (I.loadImage(FileName) != 0) {
		cout << endl << "File Error: Image Not Loaded" << endl << endl;
	}
	else {
		cout << endl << " Image has Been Loaded" << endl << endl;
	}
}

void SAVE(Image &I) {
	if (!I.isAvailable()) {
		cout << endl << "No Image Available" << endl << endl;
		return;
	}
	string FileName = "";
	cout << "Specify Image File Name: ";
	do {
		getline(cin, FileName);
	} while (FileName.length() < 1);

	int Result = I.saveImage(FileName);

	if (2 == Result) {
		cout << endl << "File Error: Image Not Saved" << endl << endl;
	}
	else {
		cout << endl << " Image Saved Successfully" << endl << endl;
	}
	return;
}

void Rotate(Image &I, Image &Result) {
	// Ask user to specify starting and ending rotation angle and a step size
	// compute and save in different files all rotated versions of the source image I.
	int CX = 0, CY = 0;
	double Angle;
    cout << "Please enter the angle: \n";
    cin >> Angle;
    cout << "\n";
	cout << "Please enter the center (x,y) , respectively : \n";
	cin >> CX >> CY;
	cout << endl;
	I.RotateImageAboutPoint(Result,Angle,CX,CY);
}

void Translate(Image &I, Image &Result) {
	// Ask user to specify Translation in X and Y and a step size
	// compute and save translated images in different files that are translated
	// versions of the source image I.
	int Tx,Ty;
	cout << "Enter the X and Y coordinates for translation : \n";
	cin >> Tx >> Ty;
	I.TranslateImage(Result,Tx,Ty);
}

void Flip(Image &I, Image &Result) {
	// Main Driver for computing Flipped version of image I
	// It must save the flipped version into image Result
	int choice;
	cout << "How you want the image to be flipped : \n\n1- Vertically ( Enter 0)\n2-Horizontally (Enter 1)\n";
	cin >> choice;
	I.FlipImage(Result, choice);
}

void Scale(Image &I, Image &Result) {
	// Main Driver for computing Flipped version of image I
	// It must save the flipped version into image Result
	double SX , SY;
	cout << "Please specify the coordinates (X and Y) : \n";
	cin >> SX >> SY;
	I.ScaleImage(Result,SX , SY);
}

void Transform(Image &I, Image &Result) {
    cout << "Please enter Transform array elements : ";
    double A[2][2];
    for(int i = 0 ; i < 2 ; i++)
    {
        for(int j = 0 ; j < 2 ; j++)
        {
            cin >> A[i][j];
        }
    }
    I.LinearTransform(Result , A);
}

void Quantize(Image &I, Image &Result ) {
	// Main Driver for Quantizing an Image
	// It must compute Quantized image by quantizing image I and save result in image Result
	int k = 26;
	int RangeArray[8][2];
	for(int i = 0 ; i < 8 ; i++)
    {
        RangeArray[i][0] = i;
    }
    for(int i = 0 ; i < 8 ; i++)
    {
        RangeArray[i][1] = k;
        k+=26;
    }
	I.QuantizeImage(Result , RangeArray);
}

void MeanFilter(Image &I, Image &Result) {
	// Main driver for mean filter
	// It must apply a 3 x 3 mean filter on image I and save result in image Result
	I.MeanFilterImage(Result,3);
}

void MedianFilter(Image &I, Image &Result) {
	// Main Driver for Median  Filter
	// It must apply 3 x 3 median filter on image I and Save Result in Image Result
	I.MedianFilterImage(Result);
}

void ComputeDerivative(Image &I, Image &Result) {
	// Main Driver for Computing Derivative
	// It must Compute Derivative of image I and Save Result in Image Result
	Mask l;
	I.ComputeDerivativeImage(Result, l);
}

void Filter(Image &I, Image &Result) {
	// Main Driver for Filter
	// It must apply a filter on image I and Save Result in Image Result
	Mask F;
	cout << "\nEnter the Mask Array elements : \n";
	for(int i = 0 ; i < 3 ; i++)
    {
        for(int j = 0 ; j < 3 ; j++)
        {
            cin >> F.MaskArray[i][j];
        }
    }
	I.ApplyFilter(Result,F);
}



int main() {
	Image A, B;
	int Choice = -2;
	Menu MainMenu("Menu.txt");

	do {
		Choice = MainMenu.ShowMenuAndGetChoice();

		if (1 == Choice) {
			LOAD(A);
		}
		else if (2 == Choice) {
			SAVE(A);
		}
		else if (3 == Choice) {
			SAVE(B);
		}
		else if (4 == Choice) {
			Rotate(A, B);
		}
		else if (5 == Choice) {
			Translate(A, B);
		}
		else if (6 == Choice) {
			Scale(A, B);
		}
		else if (7 == Choice) {
			Transform(A, B);
		}
		else if (8 == Choice) {
			Flip(A, B);
		}
		else if (9 == Choice) {
			Quantize(A, B);
		}
		else if (10 == Choice) {
			MeanFilter(A, B);
		}
		else if (11 == Choice) {
			MedianFilter(A, B);
		}
		else if (12 == Choice) {
			ComputeDerivative(A, B);
		}
		else if (13 == Choice) {
			Filter(A, B);
		}

	} while (Choice != (MainMenu.GetMaxOptions() + 1));

	cout << endl << "\n\tLeaving already. \n\tPlease come back again. \n\tSee You Soon\t";
	return 0;
}

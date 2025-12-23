#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;

class FFT{
private:
    // some simple ffts for mixed radix
    void dft_radix2 (Complex result[], Complex input[], bool inverse ){
        result[0] = input[0] + input[1];
        result[1] = input[0] - input[1];
        //special scaling for inverse
        if (inverse){
            result[0] = 0.5*result[0];
            result[1] = 0.5*result[1];
        }
    }
    void dft_radix3 (Complex result[], Complex input[], bool inverse){

    }
public:
    static bool isValidLen(size_t n){
        if (n == 0){
            return false;
        }
        while (n%2 == 0) {n /=2;}
        while (n%3 == 0) {n /=3;}
        while (n%5 == 0) {n /=5;}

        if (n == 1){
            std::cout << "Длина удволитворяет условию "<< std::endl;
        }
        else{
            std::cout << "Неверная длинна данных" << std::endl;
        }
        return n == 1;
    } 

    static void BitRev(std::vector<Complex>& data, size_t n){
        for (size_t i = 1, j = 0; i < n; ++i) {
        
        size_t bit = n >> 1;
        
        while (j & bit) {
            j ^= bit;     
            bit >>= 1;    
        }
        j ^= bit;  
        
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }  
    }
    /*static std::vector<Complex> ReadFrom(const std::string& filename){
            // INPUT : x_re , x_im Надо подумать, так как данные вещественные, если 
            std::vector<std::complex<double>> data;
            std::ifstream in("input.txt");
            double real, imag;
            while(in >> real >> imag){
                data.emplace_back(real, imag);
            }

    }*/
    static void WriteInto(const std::vector<std::complex<double>>& data, const std::string& filename){

        std::ofstream out("output.txt");
        for (const auto& value : data) {
            out << value.real() << " " << value.imag() << "\n";
            }
        }        
    
    static void fft_radix2_inplace(std::vector<Complex>& input, size_t n, bool inverse){
        if (n == 1) return;

        // optimization with pre-calc of twiddle
        std::vector<Complex> twiddle(n/2);

        for(size_t i = 0; i < n/2; i++){
            int sign = -1;
            if (inverse){sign = 1;}
            double angle = sign*2*PI*i/n;
            twiddle[i] = Complex(std::cos(angle), std::sin(angle));
        }
        BitRev(input, n);
        // main for. divide into groups 2^... 
        for (size_t len = 2; len <=n; len<<=1){
            // calculate step for pre-calc twiddle
            size_t step = n/len;

            for(size_t i = 0; i < n ; i+=len)
                for (size_t j = 0; j < len/2 ; ++j){
                    size_t idx = j*step;
                    Complex u = input[i+j];
                    Complex v = twiddle[idx]*input[i+j + len/2];

                    input[i+j] = u + v;
                    input[i+j+len/2] = u - v;
                }
        }
        
        if (inverse){
            for (size_t i = 0; i < n ; i++){ input[i] /= n;}

        }
    }

         
    
};




int main()
{
    {
        size_t n4 = 4;
        std::vector<Complex> signal = {{4,0} , {0,0} , {0,0}, {0,0}};
        FFT::fft_radix2_inplace(signal, n4, true);

        std::cout << "Результат :\n";
        for (size_t i = 0; i < 4; ++i) {
            std::cout << "X[" << i << "] = " << signal[i] << "\n";
        }
        
    }

    return 0;
}
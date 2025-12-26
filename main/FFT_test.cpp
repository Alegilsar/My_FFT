#include <iostream>
#include <random>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;

class FFT{
private:

    std::vector<Complex> input;
    size_t n;
    // some short ffts for mixed radix
    
    static std::vector<Complex> dft_radix5(std::vector<Complex>& input){
        std::vector<Complex> result(5);
        double angle  = -2*PI/5;
        Complex w = std::exp(Complex(0, angle));
        Complex w2 = w*w;
        Complex w3 = w2*w;
        Complex w4 = w2*w2;

        result[0] = input[0] +  input[1] + input[2] + input[3] + input[4];
        result[1] = input[0] +  w*input[1] + w2*input[2] + w3*input[3] + w4*input[4];
        result[2] = input[0] + w2*input[1] + w4*input[2] + w*input[3] + w3*input[4];
        result[3] = input[0] + w3*input[1] + w*input[2]  + w4*input[3] + w2*input[4];
        result[4] = input[0] + w4*input[1] + w3*input[2] + w2*input[3] + w*input[4];


        return result;
    }

    static std::vector<Complex> dft_radix2 (std::vector<Complex>& input){

        std::vector<Complex> result(2);
        // butterfly-2
        result[0] = input[0] + input[1];
        result[1] = input[0] - input[1];

        return result;
    }
    
    static std::vector<Complex> dft_radix3 (std::vector<Complex>& input){

        std::vector<Complex> result(3);
        double angle = -1*2*PI/3;

        // twiddels for radix3 but
        Complex w_3 = std::exp(Complex(0,angle));
        Complex w_3_2 = w_3*w_3;

        //butterfly-3
        result[0] = input[0] + input[1] + input[2];
        result[1] = input[0] + w_3*input[1] + w_3_2*input[2];
        result[2] = input[0] + w_3_2*input[1] + w_3*input[2];

        return result;
    }
    static bool isValidLen(size_t n){
        if (n == 0){
            return false;
        }
        while (n%2 == 0) {n /=2;}
        while (n%3 == 0) {n /=3;}
        while (n%5 == 0) {n /=5;}

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

    static std::vector<size_t> GetFactor(size_t n){
        // for mixed radix need multidim. Size of each dim will get for each radixs
        std::vector<size_t> primes = {2, 3, 5};
        std::vector<size_t> factors;
        size_t temp = n;
        
        for(size_t p:primes){
            while(temp%p == 0){
                factors.push_back(p);
                temp/=p;
            }
        }
        return factors;

    } 
    static void digit_Rev(std::vector<Complex>& input, size_t n, std::vector<size_t>& factor ){
        std::vector<size_t> reversed(n);
        for (size_t i = 0; i < n ; ++i){
            size_t rev = 0;
            size_t temp = i;
            for(int f = factor.size() - 1; f >= 0 ; --f){
                size_t radix = factor[f];
                rev = radix*rev + (temp%radix);
                temp/=radix;
            }
            reversed[i] = rev;
    
        }

        for (size_t i = 0; i < n; ++i){
            size_t j = reversed[i];
            while (j < i)
            {
                j = reversed[j];
            }
            if(j > i){
                std::swap(input[i],input[j]);
            }
        }
    
        
    }
    
    
public:
    
    
    static std::vector<Complex> fft_radix2_inplace(std::vector<Complex>& input, size_t n){
        // it was inplaced but for comfort use out of place
        std::vector<Complex> res = input;
    
        // optimization with pre-calc of twiddle
        std::vector<Complex> twiddle(n/2);

        for(size_t i = 0; i < n/2; i++){
            double angle = -1*2*PI*i/n;
            twiddle[i] = Complex(std::cos(angle), std::sin(angle));
        }
        BitRev(res, n);
        // main for. divide into groups 2^... 
        for (size_t len = 2; len <=n; len<<=1){
            // calculate step for pre-calc twiddle
            size_t step = n/len;

            for(size_t i = 0; i < n ; i+=len)
                for (size_t j = 0; j < len/2 ; ++j){
                    size_t idx = j*step;
                    Complex u = res[i+j];
                    Complex v = twiddle[idx]*res[i+j + len/2];

                    res[i+j] = u + v;
                    res[i+j+len/2] = u - v;
                }
        }
        return res;

    }
    static std::vector<Complex> fft_radix3_rec(std::vector<Complex>& input, size_t n){
        if (n==1) return input;

        // this part will be out off part because I want
        
        std::vector<Complex> X0(n/3), X1(n/3) , X2(n/3);

        for (size_t i = 0 ; i < n/3 ; i++){
            X0[i] = input[3*i];
            X1[i] = input[3*i + 1];
            X2[i] = input[3*i + 2];
        }
        //recursive for each group
        auto x0 = fft_radix3_rec(X0, n/3);
        auto x1 = fft_radix3_rec(X1, n/3);
        auto x2 = fft_radix3_rec(X2, n/3);

        std::vector<Complex> result(n);

        Complex w = std::exp(Complex(0, -2*PI/n));
        for (size_t k = 0; k < n/3 ; ++k){
            Complex w_k = std::pow(w,k);
            Complex w_k2 = w_k*w_k;

            Complex x1_twid = w_k *  x1[k];
            Complex x2_twid = w_k2 *  x2[k];

            Complex w3 = std::exp(Complex(0, -2*PI/3));
            Complex w3_2 = w3*w3;


            result[k] = x0[k] + x1_twid + x2_twid;
            result[k + n/3] = x0[k] + x1_twid*w3 + x2_twid*w3_2;
            result[k + 2*n/3] = x0[k] + x1_twid*w3_2 + x2_twid*w3;
        }
        return result;
    }
    static std::vector<Complex> mixed_radix(std::vector<Complex>& input, size_t n){
        // will be out of place, because problem with order 
        auto factors = GetFactor(n);

        std::vector<Complex> temp(n);

        std::vector<Complex>* current = &input;
        std::vector<Complex>* next = &temp;
        
        for (size_t stage = 0 ; stage < factors.size(); ++stage){
            size_t radix = factors[stage];
            size_t rem = current -> size()/ radix;

            for (size_t i = 0; i < rem ; ++i){
                std::vector<Complex> block(radix);
                for (size_t r = 0 ; r < radix; ++r ){
                    block[r] = (*current)[i + r*rem]; 
                }

                std::vector<Complex> transform(radix);
                if (radix == 2){
                    transform = dft_radix2(block);
                }
                else if(radix == 3){
                    transform = dft_radix3(block);   
                }
                else if(radix == 5){
                    transform = dft_radix5(block);
                }
                // count twiddles between stages
                for(size_t k = 0; k < radix; ++k){
                    double angle = -2*PI*i*k/n;
                    Complex twiddle = std::exp(Complex(0, angle));

                    size_t ind = i*radix + k;
                    (*next)[ind] = transform[k]*twiddle;    
                }

            }

            std::swap(current, next);
            
        }
        std::vector<Complex> result = *current;

        digit_Rev(result, n, factors);
        return result;
    }
    
    static std::vector<Complex> fft(  std::vector<Complex>& input , size_t n){
        size_t temp = n;
        std::vector<Complex> out;
        while (temp%2 == 0){
            temp/=2;
        }
        if (temp == 1){
            out = fft_radix2_inplace(input, n);
            return out ;
        }
        else{temp = n;}

        while (temp%3 == 0){
            temp/=3;
        }
        if(temp == 1){
            out = fft_radix3_rec(input, n);
            return out;
        }else{
            out = mixed_radix(input,n);
            return out;
        }
    }
    static std::vector<Complex> transform (std::vector<Complex>& input, size_t n, bool inverse ){

        bool legal =  isValidLen(n);
        if (legal){
            std::vector<Complex> out;
            if (!inverse){
                out = fft(input, n );
                return out;
            }
            else{ 
                for (size_t t = 0; t<n ; ++t){
                    input[t] = conj(input[t]);
                }
                std::vector<Complex> out = fft(input, n);
                Complex temp;
                for (size_t t = 0; t< n ; ++t){
                    temp = conj(out[t]);
                    temp /= n;
                    out[t] = temp;
                }
                return out;
            }
        }
        else{
            std::cout << "===== Invalid length =====" << std::endl;
            std::cout << "  Input array returned  " << std::endl;
            return input;
        }
    }
    
    static std::vector<Complex> generate_complex_vector(
        size_t n, 
        double real_min = -100.0, 
        double real_max = 100.0,
        double imag_min = -100.0, 
        double imag_max = 100.0) {
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> real_dist(real_min, real_max);
        std::uniform_real_distribution<> imag_dist(imag_min, imag_max);
        
        std::vector<Complex> result(n);

        
        for (size_t i = 0; i < n; ++i) {
            result[i] = (real_dist(gen), imag_dist(gen));
        }
        
        return result;

    
    } 
    static void get_error(std::vector<Complex>& input , std::vector<Complex>& out ){
        size_t N = input.size();
        Complex error;
        for (size_t k = 0; k < N; ++k){
            error  = input[k] - out[k];
            std::cout <<"Error in "<< k << " symb "<< error << std::endl;
        }
    }
};




int main()
{
    
   

    


    {
        // mixed version don't work correctly( problems with dig_rev and twiddels between stages), I'm sorry
         std::cout << "=== Тест Radix-MIxed ===\n";
        size_t n_mixed = 12;
        std::vector<Complex> test(n_mixed);
        for (size_t k = 0; k < n_mixed; ++k){
            test[k] = {1,0};
        }        
        
        std::vector<Complex> out = FFT::mixed_radix(test, n_mixed);
        std::cout<<"Полученный спектр"<< std::endl;
        for(size_t k = 0; k < n_mixed; k++){
            std::cout<< "X[" << k << "] = " << out[k] << "\n";
        }
    }   

    {
        std::cout<< "======== complited  func check ========" << "\n";
        size_t n4 = 9;
        std::vector<Complex> test = FFT::generate_complex_vector(n4) ;

        
        std::vector<Complex> out;
        out = FFT::transform(test, n4, false);

        std::vector<Complex> rev = FFT::transform(out, n4, true);

        FFT::get_error(test, rev);

    }
    return 0;


}
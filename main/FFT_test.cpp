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

    size_t n;


    // some short ffts for mixed radix
    
    std::vector<Complex> dft_radix5(std::vector<Complex>& input){
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

    std::vector<Complex> dft_radix2 (std::vector<Complex>& input){

        std::vector<Complex> result(2);
        // butterfly-2
        result[0] = input[0] + input[1];
        result[1] = input[0] - input[1];

        return result;
    }
    
    std::vector<Complex> dft_radix3 (std::vector<Complex>& input){

        std::vector<Complex> result(3);
        double angle = -2*PI/3;

        // twiddels for radix3 but
        Complex w_3 = std::exp(Complex(0,angle));
        Complex w_3_2 = w_3*w_3;

        //butterfly-3
        result[0] = input[0] + input[1] + input[2];
        result[1] = input[0] + w_3*input[1] + w_3_2*input[2];
        result[2] = input[0] + w_3_2*input[1] + w_3*input[2];

        return result;
    }
    
    // check the length
    bool isValidLen(size_t n){
        if (n == 0){
            return false;
        }
        while (n%2 == 0) {n /=2;}
        while (n%3 == 0) {n /=3;}
        while (n%5 == 0) {n /=5;}

        return n == 1;
    } 
    
    std::vector<size_t> GetFactor(size_t n){
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
    
    //Peermutation of vector components
    void digit_revers(std::vector<Complex>& input, size_t n, std::vector<size_t> factor ){
        std::vector<Complex> output(n);
        for (size_t i = 0; i < n ; ++i){
            size_t rev = 0;
            size_t temp = i;
            for(int f = factor.size() - 1 ; f >= 0 ; --f){
                size_t radix = factor[f];
                rev = radix*rev + (temp%radix);
                temp/=radix;
            }
        
            output[rev] = input[i];
        }
        input = std::move(output);
        
    }
 
    void BitRev(std::vector<Complex>& data, size_t n){
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
public:
    
    std::vector<Complex> fft_radix2(std::vector<Complex>& input, size_t n){
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
    std::vector<Complex> fft_radix3_rec(std::vector<Complex>& input, size_t n){
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
 
 
    void process_stage(std::vector<Complex>& data, size_t stage, size_t stage_size, std::vector<size_t> factors, size_t N){
        size_t radix = factors[stage];
        size_t groups = N / (radix * stage_size);
    
        std::vector<Complex> block(radix);
        
        for (size_t group = 0; group < groups; ++group) {
            for (size_t offset = 0; offset < stage_size; ++offset) {


                // Сonstruct block + apply twiddles 
                size_t baseIdx = group * radix * stage_size + offset;
                for (size_t r = 0; r < radix; ++r) {
                    size_t idx = baseIdx + r * stage_size;
                    double angle = -2.0 * PI * r * offset *groups /(N);

                    block[r] = data[idx]*std::exp(Complex(0, angle));
                }
                //Calculate butterfly of block 
                std::vector<Complex> dftResult(radix);
                switch (radix) {
                    case 2:
                        dftResult = dft_radix2(block);
                        break;
                    case 3:
                        dftResult = dft_radix3(block);
                        break;
                    case 5:
                        dftResult = dft_radix5(block);
                        break;

                }
            
                // rewrite input with dft_result
                for (size_t r = 0; r < radix; ++r) {
                    size_t idx = baseIdx + r * stage_size;
                    data[idx] = dftResult[r];
                }
            }
        }

    }

    std::vector<Complex> mixed_radix (std::vector<Complex>& input , size_t n ){
     
        auto factors = GetFactor(n);

        digit_revers (input, n, factors);
        size_t stage_size = 1;
        for (size_t stage = 0; stage < factors.size(); ++stage){
            process_stage(input , stage, stage_size, factors, n);
            stage_size *= factors[stage];

        }
        return input;
    }

    std::vector<Complex> fft(  std::vector<Complex>& input , size_t n){
        size_t temp = n;
        std::vector<Complex> out;
        while (temp%2 == 0){
            temp/=2;
        }
        if (temp == 1){
            out = fft_radix2(input, n);
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
    std::vector<Complex> transform (std::vector<Complex>& input, size_t n, bool inverse ){

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
    
    
};

std::vector<Complex> generate_complex_vector(
    size_t len, 
    double real_min = -100.0, 
    double real_max = 100.0,
    double imag_min = -100.0, 
    double imag_max = 100.0) {
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> real_dist(real_min, real_max);
    std::uniform_real_distribution<> imag_dist(imag_min, imag_max);
    
    std::vector<Complex> result(len);

    
    for (size_t i = 0; i < len; ++i) {
        result[i] = {real_dist(gen), imag_dist(gen)};
    }
    
    return result;


} 
void get_error(std::vector<Complex>& input , std::vector<Complex>& out ){
        size_t N = input.size();
        Complex error;
        for (size_t k = 0; k < N; ++k){
            error  = input[k] - out[k];
            std::cout <<"Error in "<< k << " symb "<< error << std::endl;
        }
}

std::vector<Complex> generate_multi_har_signal(std::vector<int>& harmonics, size_t length)
{
    
    std::cout<< "Используемые гармоники :\n";
    for (size_t i = 0 ; i < harmonics.size() ; ++i){
         std::cout<< harmonics[i]<< "\n" ;
    }
  
    // generate multi-harmonic signal 
    std::vector<Complex> signal(length);
    for(size_t k = 0 ; k < length ; ++k){
        for(size_t i = 0 ; i < harmonics.size(); ++i){
            int har = harmonics[i];
            double angle = 2*PI*har*k/length;
            signal[k] += std::exp(Complex(0,angle));
        }        
    }

    return signal;
}  

int main()
{
    

    {
        std::cout << "==== Test n-harmonics for radix-2 ===="<< std::endl;
        size_t n2 = 8;
        
        std::vector<int> components = {3, 7};      
        std::vector<Complex> sig = generate_multi_har_signal(components, n2);
        FFT F2;
        
        std::vector<Complex> out_2 = F2.fft_radix2(sig, n2);

        std::cout<<"Полученный спектр для radix-2"<< std::endl;
        for(size_t k = 0; k < n2; k++){
            std::cout<< "X[" << k << "] = " << out_2[k] << "\n";
        }

    }

    {
        std::cout << "==== Test n-harmonics for radix-3 ===="<< std::endl;
        size_t n3 = 9;
        
        std::vector<int> components = {3, 7};      
        std::vector<Complex> sig = generate_multi_har_signal(components, n3);
        FFT F3;
        
        std::vector<Complex> out = F3.fft_radix3_rec(sig, n3);

        std::cout<<"Полученный спектр для radix-3"<< std::endl;
        for(size_t k = 0; k < n3; k++){
            std::cout<< "X[" << k << "] = " << out[k] << "\n";
        }

    }

    {
        std::cout << "==== Test n-harmonics for radix-5 ===="<< std::endl;
        size_t n5 = 25;
        
        std::vector<int> components = {3, 7};      
        std::vector<Complex> sig = generate_multi_har_signal(components, n5);
        FFT F5;
        
        std::vector<Complex> out = F5.mixed_radix(sig, n5);

        std::cout<<"Полученный спектр для radix-5"<< std::endl;
        for(size_t k = 0; k < n5; k++){
            std::cout<< "X[" << k << "] = " << out[k] << "\n";
        }

    }

    {
        
       std::cout << "==== Test n-harmonics for mixed - radix ===="<< std::endl;
        size_t n_mixed = 24;
        
        std::vector<int> components = {3, 7};      
        std::vector<Complex> sig = generate_multi_har_signal(components, n_mixed);
        FFT F3;
        
        std::vector<Complex> out = F3.mixed_radix(sig, n_mixed);

        std::cout<<"Полученный спектр для mixed radix"<< std::endl;
        for(size_t k = 0; k < n_mixed; k++){
            std::cout<< "X[" << k << "] = " << out[k] << "\n";
        }

    }


    {
        std::cout<< "======== complited  func check ========" << "\n";
        size_t n4 = 12;
        FFT mixed_test;
        std::vector<Complex> test = generate_complex_vector(n4) ;
        std::vector<Complex> test_copy = test;

        
        std::vector<Complex> spectrum =  mixed_test.transform(test_copy, n4, false);
        std::vector<Complex> time_output =  mixed_test.transform(spectrum, n4, true);


        get_error(test, time_output );

    }
    return 0;


}
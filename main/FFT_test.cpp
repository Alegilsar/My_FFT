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
    std::vector<std::vector<Complex>> twiddleTables;
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
    
    
    
public:
    
    static void digit_Rev(std::vector<Complex>& input, size_t n, std::vector<size_t> factor ){
        std::vector<Complex> output(n);
        for (size_t i = 0; i < n ; ++i){
            size_t rev = 0;
            size_t temp = i;
            for(int f = 0; f < factor.size() ; ++f){
                size_t radix = factor[f];
                rev = radix*rev + (temp%radix);
                temp/=radix;
            }
        
            output[rev] = input[i];
        }
        input = std::move(output);
        
    }
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
                    std::cout << twiddle[idx] << std::endl;
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
 
    static void process_stage(std::vector<Complex>& data, size_t stage, size_t stage_size, std::vector<size_t> factors, size_t N, std::vector<std::vector<Complex>>& twiddleTables){
        size_t radix = factors[stage];
        size_t groups = N / (radix * stage_size);
        std::cout<< " =======Cтадияы == " << stage << std::endl;
        const auto& twiddles = twiddleTables[stage];
        size_t twiddleIdx = 0;
        
    
        std::vector<Complex> block(radix);
        
        for (size_t group = 0; group < groups; ++group) {
            for (size_t offset = 0; offset < stage_size; ++offset) {


                // 1. Собираем блок из данных

                size_t baseIdx = group * radix * stage_size + offset;
                for (size_t r = 0; r < radix; ++r) {
                    size_t idx = baseIdx + r * stage_size;
                    block[r] = data[idx];
                }
                for (size_t r = 0; r < radix; ++r) {
                    size_t k = offset;

                    double angle = -2.0 * PI * r * k *groups /(N);
                    std::cout <<"Индекс"<< std::endl;
                    std::cout << r*k % N << std::endl;
                    std::cout << std::exp(Complex(0, angle)) << std::endl;
                    block[r] *= std::exp(Complex(0, angle));
                }
                std::vector<Complex> dftResult;
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
                    default:
                        throw std::runtime_error("Unsupported radix");
                }
                
                // 3. Применяем twiddle factors (кроме первого элемента)
               
                
                // 4. Записываем результат обратно
                for (size_t r = 0; r < radix; ++r) {
                    size_t idx = baseIdx + r * stage_size;
                    data[idx] = dftResult[r];
                }
            }
        }

    }

    static std::vector<Complex> mixed_radix (std::vector<Complex>& input , size_t n ){
        auto factors = GetFactor(n);
        std::vector<std::vector<Complex>> twiddleTables;
        twiddleTables.clear();
        size_t stageSize = 1;
        
        for (size_t stage = 0; stage < factors.size(); ++stage) {
            
            size_t radix = factors[stage];
            size_t groups = n / (radix * stageSize);
            std::vector<Complex> stageTwiddles;
            
            if (stage  == 0){
                for (size_t group = 0; group < groups; ++group)
                 stageTwiddles.push_back(std::exp(Complex(0, 0)));
            }
            else{
            // Для каждого элемента внутри блока (кроме первого)
            for (size_t group = 0; group < groups; ++group) {
                for (size_t k = 0; k < stageSize; ++k) {
                    double baseAngle = -2.0 * PI * (k )/ (radix*stageSize);
                    
                    for (size_t r = 1; r < radix; ++r) {
                        double angle = baseAngle * r;
                        stageTwiddles.push_back(std::exp(Complex(0, angle)));
                    }
                }
            }
            }
            
            twiddleTables.push_back(stageTwiddles);
            stageSize *= radix;
        }
        
        digit_Rev(input, n, {5,2});
        size_t stage_size = 1;
        for (size_t stage = 0; stage < factors.size(); ++stage){
            process_stage(input , stage, stage_size, factors, n, twiddleTables);
            stage_size *= factors[stage];

        }
        return input;
    }

    /*static void precomputeTwiddles() {
        twiddleTables.clear();
        size_t stageSize = 1;
        
        for (size_t radix : factors) {
            size_t groups = N / (radix * stageSize);
            std::vector<Complex> stageTwiddles;
            
            // Twiddle factors для этого этапа
            // Для каждого элемента внутри блока (кроме первого)
            for (size_t group = 0; group < groups; ++group) {
                for (size_t k = 0; k < stageSize; ++k) {
                    double baseAngle = -2.0 * PI * (group * stageSize + k) / N;
                    
                    for (size_t r = 1; r < radix; ++r) {
                        double angle = baseAngle * r;
                        stageTwiddles.push_back(std::exp(Complex(0, angle)));
                    }
                }
            }
            
            twiddleTables.push_back(stageTwiddles);
            stageSize *= radix;
        }
    }*/
    /*   static std::vector<Complex> mixed_radix(std::vector<Complex>& input, size_t n){
        // will be out of place, because problem with order 
        auto factors = GetFactor(n);

        std::vector<Complex> result  = input;
        digit_Rev(result, n, factors);

        
        for (size_t stage = 0 ; stage < factors.size(); ++stage){
            size_t radix = factors[stage];
            size_t rem = current -> size()/ radix;

            for (size_t i = 0; i < rem ; ++i){
                std::vector<Complex> block(radix);
                for (size_t r = 0 ; r < radix; ++r ){
                    double angle = -2*PI*i*r/n;
                    Complex twiddle = std::exp(Complex(0, angle));

                    block[r] = (*current)[i + r*rem]*twiddle; 
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
                    
                    size_t ind = i*radix + k;
                    (*next)[ind] = transform[k];    
                }

                

            }

            std::swap(current, next);
            
        }
        std::vector<Complex> result = *current;

        //digit_Rev(result, n, factors);
        return result;
    }*/
    
    static void Stockham_stage (std::vector<Complex>& input, std::vector<Complex>& output, size_t N, size_t radix, size_t stage_size){

        size_t groups = N/(radix*stage_size);

        for (size_t group = 0; group < groups; ++group) {
            for (size_t offset = 0; offset < stage_size; ++offset) {
                std::vector<Complex> temp(radix);

                for (size_t r = 0; r < radix; ++r) {
                    size_t idx = group*radix*stage_size + r * stage_size + offset;
                    size_t g_k = group*stage_size+ offset;
                    double angle = -2.0*PI*r*offset*groups/(N);
                    temp[r] = input[idx]*std::exp(Complex(0, angle));
                }
                
                std::vector<Complex> dftResult;
                switch (radix) {
                    case 2: dftResult = dft_radix2(temp); break;
                    case 3: dftResult = dft_radix3(temp); break;
                    case 5: dftResult = dft_radix5(temp); break;
                    default:
                        throw std::runtime_error("Unsupported radix");
                }
                
                // twiddle + place
                for (size_t r = 0; r < radix; ++r){
                    size_t out_ind = group*radix*stage_size + offset*radix + r;
                    
                    output[out_ind] = dftResult[r];
                }
               
            }
        }

    }
    static std::vector<Complex> stockham_mixed_radix_fft(const std::vector<Complex>& input) {
        
        size_t N = input.size();
        std::vector<size_t>factors = {2,2,2};

        std::vector<Complex> buf1 = input;
        std::vector<Complex> buf2(N);

        size_t stageSize = 1;
        bool bufer_switch = true;

        for (size_t stage = 0; stage < factors.size(); ++stage) {
            size_t radix = factors[stage];

            if (bufer_switch)
                Stockham_stage(buf1, buf2, N, radix, stageSize);
            else
                Stockham_stage(buf2, buf1, N, radix, stageSize);

            stageSize *= radix;
            bufer_switch = !bufer_switch;
        }

        return bufer_switch ? buf1 : buf2;
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
            //out = mixed_radix(input,n);
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

void fft_stage_order(std::vector<Complex>& input , std::vector<size_t> factor, size_t N ){
    

}


int main()
{
    
    /*{
        std::vector<size_t> factors = {2,2,3};
        std::vector<Complex> sig = { {0,0}, {1,0} , {2,0} , {3,0} , {4,0} , {5,0}, {6, 0}, {7,0}, {8,0}, {9,0},{10,0}, {11,0}};
        FFT::digit_Rev(sig, 12, factors);
        size_t N = sig.size();
        std::vector<Complex> copy = sig;
        size_t stage_size = 1;
        for (auto f : factors){
            size_t radix = f;
            size_t groups = N / (radix*stage_size);
            for (size_t group = 0; group < groups ; ++group ){
                for (size_t offset = 0; offset < stage_size ; ++offset  ){
                        size_t base_ind = group * radix * stage_size + offset;
                        std::vector<Complex> block(radix);
                        for (size_t r = 0; r < radix ; ++r){
                            size_t idx = base_ind + r * stage_size;
                            block[r] = sig[idx];
                        }  

                        for (size_t r = 0; r < radix; ++r) {
                            size_t idx = base_ind*stage_size + r;
                            copy[idx] = block[r];
                        }
                }
            }

            stage_size *= radix;

        }

    }*/

    {
        std::cout << "Comparison of two methods radix_2 and mixed-radix for len = 2^k"<< std::endl;
        size_t n2 = 8;
        std::vector<Complex> sig(n2);
        for (size_t k = 0; k < n2; ++k){
            sig[k] = std::exp(Complex(0, 2*PI*3*k/n2)) + std::exp(Complex(0, 2*PI*6*k/n2));
        }       
        
        std::vector<Complex> out_2 = FFT::fft_radix2_inplace(sig, n2);

        std::cout<<"Полученный спектр для radix-2"<< std::endl;
        for(size_t k = 0; k < n2; k++){
            std::cout<< "X[" << k << "] = " << out_2[k] << "\n";
        }

    }
    


 {
        
        std::cout << "=== Тест Radix-3 ===\n";
        size_t n_mixed = 10;
        std::vector<Complex> test(n_mixed);
        for (size_t k = 0; k < n_mixed; ++k){
           test[k] =  std::exp(Complex(0, 2*PI*3*k/n_mixed));
        }
        /*test[0] =Complex(1.0,0.0); 
        for (size_t k = 1; k < n_mixed; ++k){

            test[k] = Complex(0,0);
        }   */     
        
        std::vector<Complex> out = FFT::mixed_radix(test, n_mixed);
        std::cout<<"Полученный спектр"<< std::endl;
        for(size_t k = 0; k < n_mixed; k++){
            std::cout<< "X[" << k << "] = " << out[k] << "\n";
        }
    }
 /*
    {
        // mixed version don't work correctly( problems with dig_rev and twiddels between stages), I'm sorry
         std::cout << "=== Тест Radix-Mixed ===\n";
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
        size_t n4 = 12;
        std::vector<Complex> test = FFT::generate_complex_vector(n4) ;

        
        std::vector<Complex> out;
        out = FFT::transform(test, n4, false);

        std::vector<Complex> rev = FFT::transform(out, n4, true);

        FFT::get_error(test, rev);

    }*/
    return 0;


}
#include <iostream>
#include <cstring>
#include <cstddef>
#include <gmp.h>

int invalid_params() {
    std::cout << "Invalid params" << std::endl;
    return 1;
}

void extended_gcd(mpz_t arg_res_r, mpz_t arg_res_t, mpz_t arg_a, mpz_t arg_b) {
    int res;
    mpz_t r, r_prev, q, t, t_prev, tmp;
    mpz_inits(q, t_prev, tmp, NULL);
    mpz_init_set(r, arg_b);
    mpz_init_set(r_prev, arg_a);
    mpz_init_set_ui(t, 1);

    while(mpz_sgn(r) != 0) {
        mpz_tdiv_qr(q, tmp, r_prev, r);
        mpz_set(r_prev, r);
        mpz_set(r, tmp);

        mpz_submul(t_prev, q, t);
        mpz_set(tmp, t);
        mpz_set(t, t_prev);
        mpz_set(t_prev, tmp);
    }

    mpz_set(arg_res_r, r_prev);
    mpz_mod(arg_res_t, t_prev, arg_a);
}

bool coprime(mpz_t a, mpz_t b) {
    mpz_t n, t;
    mpz_inits(n, t, NULL);
    extended_gcd(n, t, a, b);
    return mpz_cmp_ui(n, 1) == 0;
}

int jacobi(mpz_t arg_a, mpz_t arg_n) {
    mpz_t none, a, n;
    int j = 1;
    mpz_init(none);
    mpz_init_set(a, arg_a);
    mpz_init_set(n, arg_n);

    while(mpz_sgn(a) != 0) {
        while(mpz_even_p(a)) {
            mpz_div_ui(a, a, 2);
            int res = mpz_mod_ui(none, n, 8);

            if(res == 3 or res == 5)
                j *= -1;
        }
        mpz_swap(a, n);

        if(mpz_mod_ui(none, a, 4) == 3 and mpz_mod_ui(none, n, 4) == 3)
            j *= -1;

        mpz_mod(a, a, n);
    }
    if(mpz_cmp_ui(n, 1) == 0)
        return j;

    return 0;
}

bool prime(gmp_randstate_t rstate, mpz_t n, int k) {
    mpz_t n_2, a, x, power, res;
    mpz_init(n_2);
    mpz_init(a);
    mpz_init(x);
    mpz_init(power);
    mpz_init(res);

    if(mpz_cmp_ui(n, 1) == 0 or mpz_cmp_ui(n, 3) == 0)
        return true;
    
    if(mpz_even_p(n))
        return false;

    mpz_sub_ui(n_2, n, 2);
    mpz_sub_ui(power, n, 1);
    mpz_tdiv_q_ui(power, power, 2);

    for(int i = 0; i < k; ++i) {
        mpz_urandomm(a, rstate, n_2);
        mpz_add_ui(a, a, 2);
        mpz_set_si(x, jacobi(a, n));
        mpz_mod(x, x, n);
        mpz_powm(res, a, power, n);

        if(mpz_sgn(x) == 0 or mpz_cmp(res, x) != 0)
            return false;
    }

    return true;
}

int generate(int argc, char **argv) {
    char *endptr;
    mpz_t p, q, d, e, n, m, phi, p_1, q_1, tmp;
    gmp_randstate_t rstate; 
    int bits = strtol(argv[2], &endptr, 10);

    mpz_inits(p, q, d, e, n, m, phi, p_1, q_1, tmp, NULL);
    gmp_randinit_default(rstate);

    do {
        mpz_urandomb(p, rstate, bits/2); // P
        mpz_setbit(p, bits/2 - 1);
    } while(!prime(rstate, p, bits/2));

    do {
        mpz_urandomb(q, rstate, bits/2); // Q
        mpz_setbit(q, bits/2 - 1);
    } while(!prime(rstate, q, bits/2));

    mpz_mul(n, p, q); // N = P * Q
    mpz_sub_ui(p_1, p, 1); // P_1 = P - 1
    mpz_sub_ui(q_1, q, 1); // Q_1 = Q - 1
    mpz_mul(phi, p_1, q_1); // PHI = (P - 1) * (Q - 1)
    mpz_sub_ui(tmp, phi, 1); // tmp = PHI - 1

    do {
        mpz_urandomm(e, rstate, tmp); // E = random(0, PHI - 2)
        mpz_add_ui(e, e, 1); // E = random(1, PHI - 1)
    } while(!coprime(e, phi));

    extended_gcd(tmp, d, phi, e);
    gmp_printf("%#Zx %#Zx %#Zx %#Zx %#Zx\n", p, q, n, e, d);

    return 0;
}

int crypt(mpz_t key, mpz_t modulus, mpz_t message) {
    mpz_t res;
    mpz_init(res);

    mpz_powm(res, message, key, modulus);
    gmp_printf("%#Zx\n", res);

    return 0;
}

int main(int argc, char **argv) {
    int (f)(mpz_t, mpz_t, mpz_t);
    mpz_t key, modulus, message;
    char option;

    if(argc < 3 or strlen(argv[1]) != 2 or argv[1][0] != '-') {
        return invalid_params();
    }

    option = argv[1][1]; // g/e/d

    switch(option) {
        case 'e':;
        case 'd': break;
        case 'g': return generate(argc, argv);
        default: return invalid_params();
    }

    if(argc != 5) 
        return invalid_params();

    mpz_inits(key, modulus, message, NULL);
    mpz_set_str(key, argv[2] + 2, 16);
    mpz_set_str(modulus, argv[3] + 2, 16);
    mpz_set_str(message, argv[4] + 2, 16);

    return crypt(key, modulus, message);
}
/*
* Created on 2026.05.23
* Copyright (c) Youcef Lemsafer
*/


print_hex(x) = {
    if(x == 0,
        print1("0x0"),
        print1("0x");
        [u | u <- digits(x, 16), printf("%x", u)];
    );
}

print_u128(x) = {
    my(x_hi = x \ 2^64);
    if(x_hi > 0,
        print1("(__uint128_t{");
        print_hex(x_hi);
        print1("ull");
        print1("} << 64) | ");
    );
    print_hex(x % 2^64);
    print1("ull");
}

gen_mu64_testcases() = {
    my(cases=[3, 5, 7, 9, 15, 17, 19, 121, 257, 521,
          991, 65535, 65537, 991*997, 2^32-5, (2^32-5)*(2^32-5),
          (2^32-1)*(2^32-1), 2^64 - 59]);
    for(i=1, #cases,
        print1("{ ", cases[i], "ull, ");
        print_hex(2^64 \ cases[i]);
        print("ull },");
    )
}

gen_mu128_testcases() = {
    my(cases=[3, 5, 7, 9, 15, 17, 19, 121, 255, 257, 521,
                991, 65535, 65537, 991*997, 2^32-5, 2^32-1,
                2^32+1, (2^32-5)*(2^32-5), (2^32-1)*(2^32-1),
                2^64 - 59, 2^64-1]);
    for(i = 1, #cases,
        print1("{ ", cases[i], "ull, ");
        print_u128(2^128 \ cases[i]);
        print(" },");
    )
}


gen_mul128hi_testcases() = {
    my(cases=[
        [0, 1],
        [1, 0],
        [2, 4],
        [2, 9],
        [4, 2],
        [9, 2],
        [4095, 2^24],
        [2^43-1, 2^43+1],
        [2^64-1, 2^64-1],
        [2^64-1, 2^64],
        [2^64, 2^64-1],
        [2^64, 2^64],
        [2^65-2^33+1, 2^65-1],
        [2^53, 2^75],
        [2^77, 2^51],
        [2^99+2^94+9876543210123456789, 2^101 + 2^37 - 2^29 + 1],
        [0x0eb7fb06849c5843071dcbfb976e9a33, 0x5cf3a89821482e95e382ea861ae682b2],
        [2^125-1, 2^67],
        [0x984f9178540196a9f64fdb69582b2964, 0xbc529db2d7a28b58da7b77de9bb5098d],
        [0xabb62948a7fb7b3f8e6cb2eb530a1cda, 0x9b6ad2d7834f3926b156d7c99f3a58d7],
        [0x668eea7b862b221f82813d5aa020b600, 0xf99246d69d14700e116bf17a60eed1c9],
        [0xf0370ab5d6a70b0492d71d59775f32bb, 0x597c74ada5452af83c4b13226d6de952],
        [0x613f90a3232a35244391fb0150ddc03f, 0xc7852de57fde44ec98cd1766b758b5aa],
        [0xfdafdafdafdafdafdafdafdafdafdafd, 0xafdafdafdafdafdafdafdafdafdafdaf],
        [0xfdfdfdfdfdfdfdfdfdfdfdfdfdfdfdfd, 0xfdfdfdfdfdfdfdfdfdfdfdfdfdfdfdfd],
        [2^128-159, 2^128-1],
        [2^128-1, 2^128-159],
        [2^128-159, 2^128-159],
        [2^128-1, 2^128-1]
        ]);
    for(i = 1, #cases,
        print1("{ ");
        my(a = cases[i][1]);
        my(b = cases[i][2]);
        print_u128(a);
        print1(", ");
        print_u128(b);
        print1(",");
        my(ab_hi = (a * b) \ 2^128);
        if(ab_hi > 2^64,
            print1("\n    "),
            print1(" ")
        );
        print_u128(ab_hi);
        print(" },");
    )
}



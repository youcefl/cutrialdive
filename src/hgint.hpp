/*
* Youcef Lemsafer
* 2023.01.18
*/
#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <exception>
#include <gmp.h>

class InvalidNumberException : public std::exception
{
public:
    InvalidNumberException(std::string&& msg)
        : what_(std::move(msg))
    {}
    char const* what() const noexcept override {
        return what_.data();
    }
private:
    std::string what_;
};

class HgInt
{
public:
    ~HgInt();
    HgInt();
    HgInt(HgInt const &);
    HgInt(HgInt&&);
    HgInt & operator=(HgInt&&);
    HgInt(std::string const &);
    HgInt(char const*);
    HgInt(uint64_t);
    bool operator!() const;
    HgInt & operator+=(HgInt const& x);
    HgInt & operator+=(uint64_t x);
    HgInt & operator-=(HgInt const& x);
    HgInt & operator-=(uint64_t x);
    HgInt & operator*=(HgInt const& x);
    HgInt & operator*=(uint64_t x);
    HgInt & operator/=(HgInt const& x);
    HgInt & operator/=(uint64_t x);
    HgInt & operator<<=(uint64_t x);
    mpz_t const & get() const;
    mpz_t & get();
    mp_limb_t const * limbs() const;
    std::size_t size() const;
    std::size_t sizeInBase(int base) const;
    friend HgInt pow(HgInt x, uint64_t y);
    uint64_t mod(uint64_t x) const;
    friend void div_mod(HgInt const & number, HgInt const & divisor, HgInt& quotient, HgInt& remainder);
    friend HgInt mod(HgInt x, uint64_t y);
    friend std::string to_string(HgInt const & x);
private:
    mpz_t val_;
};

bool is_prp(HgInt number);

#define HGINT_INLINE inline

HGINT_INLINE
mpz_t const & HgInt::get() const
{
    return val_;
}

HGINT_INLINE
mpz_t & HgInt::get()
{
    return val_;
}

HGINT_INLINE
mp_limb_t const * HgInt::limbs() const
{
    return mpz_limbs_read(val_);
}

HGINT_INLINE
std::size_t HgInt::size() const
{
    return mpz_size(val_);
}

HGINT_INLINE
std::size_t HgInt::sizeInBase(int base) const
{
    return mpz_sizeinbase(val_, base);
}


HGINT_INLINE
HgInt::~HgInt()
{
    mpz_clear(val_);
}

HGINT_INLINE
HgInt::HgInt()
{
    mpz_init(val_);
}

HGINT_INLINE
HgInt::HgInt(HgInt const & other)
{
    mpz_init(val_);
    mpz_set(val_, other.val_);
}
HGINT_INLINE
HgInt::HgInt(HgInt&& other)
{
    mpz_init(val_);
    mpz_swap(val_, other.val_);
}
HGINT_INLINE
HgInt & HgInt::operator=(HgInt&& other)
{
    if (this != &other) {
        mpz_set(val_, other.val_);
    }
    return *this;
}

HGINT_INLINE
HgInt::HgInt(uint64_t v)
{
    mpz_init_set_ui(val_, v);
}

HGINT_INLINE
HgInt::HgInt(std::string const & str)
    : HgInt(str.c_str())
{
}

HGINT_INLINE
HgInt::HgInt(char const* str)
{
    if(mpz_init_set_str(val_, str, 0) != 0) {
        mpz_clear(val_);
        throw InvalidNumberException(str);
    }
}

HGINT_INLINE
bool HgInt::operator!() const {
    return mpz_sgn(val_) == 0;
}

HGINT_INLINE
HgInt & HgInt::operator+=(HgInt const& x)
{
    mpz_add(val_, val_, x.val_);
    return *this;
}
HGINT_INLINE
HgInt & HgInt::operator+=(uint64_t x)
{
    mpz_add_ui(val_, val_, x);
    return *this;
}
HGINT_INLINE
HgInt operator+(HgInt x, HgInt const & y)
{
    return x += y;
}
HGINT_INLINE
HgInt operator+(HgInt x, uint64_t y)
{
    return x += y;
}


HGINT_INLINE
HgInt & HgInt::operator-=(HgInt const& x)
{
    mpz_sub(val_, val_, x.val_);
    return *this;
}
HGINT_INLINE
HgInt & HgInt::operator-=(uint64_t x)
{
    mpz_sub_ui(val_, val_, x);
    return *this;
}
HGINT_INLINE
HgInt operator-(HgInt x, HgInt const & y)
{
    return x -= y;
}
HGINT_INLINE
HgInt operator-(HgInt x, uint64_t y)
{
    return x -= y;
}

HGINT_INLINE
HgInt & HgInt::operator*=(HgInt const& x)
{
    mpz_mul(val_, val_, x.val_);
    return *this;
}
HGINT_INLINE
HgInt & HgInt::operator*=(uint64_t x)
{
    mpz_mul_ui(val_, val_, x);
    return *this;
}
HGINT_INLINE
HgInt operator*(HgInt x, HgInt const & y)
{
    return x *= y;
}
HGINT_INLINE
HgInt operator*(HgInt x, uint64_t y)
{
    return x *= y;
}

HGINT_INLINE
HgInt & HgInt::operator/=(HgInt const& x)
{
    mpz_fdiv_q(val_, val_, x.val_);
    return *this;
}
HGINT_INLINE
HgInt & HgInt::operator/=(uint64_t x)
{
    mpz_fdiv_q_ui(val_, val_, x);
    return *this;
}
HGINT_INLINE
HgInt operator/(HgInt x, HgInt const & y)
{
    return x /= y;
}
HGINT_INLINE
HgInt operator/(HgInt x, uint64_t y)
{
    return x /= y;
}

HGINT_INLINE
HgInt & HgInt::operator<<=(uint64_t x)
{
    mpz_mul_2exp(val_, val_, x);
    return *this;    
}
HGINT_INLINE
HgInt operator<<(HgInt x, uint64_t y)
{
    return x <<= y;
}

HGINT_INLINE
HgInt pow(HgInt x, uint64_t y)
{
    mpz_pow_ui(x.val_, x.val_, y);
    return x;
}

HGINT_INLINE
uint64_t HgInt::mod(uint64_t x) const
{
    return mpz_fdiv_ui(val_,  x);
}

HGINT_INLINE
void div_mod(HgInt const& number, HgInt const& divisor, HgInt& quotient, HgInt& remainder)
{
    mpz_fdiv_qr(quotient.val_, remainder.val_, number.val_, divisor.val_);
}

HGINT_INLINE
HgInt mod(HgInt x, uint64_t y)
{
    return HgInt{x.mod(y)};
}

HGINT_INLINE
std::ostream&
operator<<(std::ostream& os, HgInt const& hg) {
    char* str = mpz_get_str(nullptr, 10, hg.get());
    os << str;
    free(str);
    return os;
}

HGINT_INLINE
std::string
to_string(HgInt const& hg) {
    char* str = mpz_get_str(nullptr, 10, hg.get());
    std::string v{str};
    free(str);
    return v;
}


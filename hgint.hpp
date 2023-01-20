/*
* Youcef Lemsafer
* 2023.01.18
*/
#pragma once

//#include <cstdint>
#include <string>
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
    HgInt & operator+=(HgInt const& x);
    HgInt & operator+=(uint64_t x);
    HgInt & operator-=(HgInt const& x);
    HgInt & operator-=(uint64_t x);
    HgInt & operator*=(HgInt const& x);
    HgInt & operator*=(uint64_t x);
    HgInt & operator/=(HgInt const& x);
    HgInt & operator/=(uint64_t x);
    mpz_t const & get() const;
    mpz_t & get();
    mp_limb_t const * limbs() const;
    std::size_t size() const;
    std::size_t sizeInBase(int base) const;
    friend HgInt pow(HgInt x, uint64_t y);
private:
    mpz_t val_;
};

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
HgInt pow(HgInt x, uint64_t y)
{
    mpz_pow_ui(x.val_, x.val_, y);
    return x;
}


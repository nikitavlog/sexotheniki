#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Simple example of encrypting, sending and decrypting a message using RSA.
Also included is how a third party without permission, would need
to use a very time consuming (brute force) approach to decrypt
the message.  This is done by finding the two prime numbers that make up
part of the public key so that you can read the encrypted message.
Implemented in python.
forked from: https://gist.github.com/JonCooperWorks/5314103
'''

import math
import random
from datetime import datetime


'''
Euclid's algorithm for determining the greatest common divisor
Use iteration to make it faster for larger integers
'''
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

'''
Euclid's extended algorithm for finding the multiplicative inverse of two numbers
'''
def multiplicative_inverse(e, phi):
    f = 0
    x1 = 0
    x2 = 1
    y1 = 1
    temp_phi = phi

    while e > 0:
        temp1 = temp_phi/e
        temp2 = temp_phi - (temp1 * e)
        temp_phi = e
        e = temp2

        x = x2- temp1 * x1
        y = f - temp1 * y1

        x2 = x1
        x1 = x
        f = y1
        y1 = y

    if temp_phi == 1:
        return f + phi


'''
Tests to see if a number is prime.
'''
def is_prime(num):
    if num in (2, 3):
        return True
    if num < 2 or num % 2 == 0:
        return False
    for n in range(3, int(num**0.5)+2, 2):
        if num % n == 0:
            return False
    return True


def generate_keypair(p, q, e=None):
    if not (is_prime(p) and is_prime(q)):
       raise ValueError('Both numbers must be prime.')
    elif p == q:  # does this matter?
        raise ValueError('p and q cannot be equal')
    # n = pq
    n = p * q

    # Phi is the totient of n
    phi = (p-1) * (q-1)

    if e is None:
        # Choose an integer e such that e and phi(n) are coprime
        # Use Euclid's Algorithm to verify that e and phi(n) are comprime
        haveCoPrime = False
        while not haveCoPrime:
            e = random.randrange(1, phi)
            g = gcd(e, phi)
            haveCoPrime = (g == 1)

    # Use Extended Euclid's Algorithm to generate the private key
    d = multiplicative_inverse(e, phi)

    # Return public and private keypair
    # Public key is (e, n) and private key is (d, n)
    return ((e, n), (d, n))


def encrypt(pk, plaintext):
    # Unpack the key into it's components
    e, n = pk
    # Convert each letter in the plaintext to numbers based on the character using a^b mod m
    #cipher = [(ord(char) ** e) % n for char in plaintext]
    cipher = [pow(ord(char), e, n) for char in plaintext]
    # Return the array of bytes
    return cipher


def decrypt(pk, ciphertext):
    # Unpack the key into its components
    d, n = pk
    # Generate the plaintext based on the ciphertext and key using a^b mod m
    #plain = [chr((char ** d) % n) for char in ciphertext]
    plain = [chr(pow(char, d, n)) for char in ciphertext]
    # Return the array of bytes as a string
    return ''.join(plain)


def primeFactorizationV1(n):
    """
    Finds the prime factors of `n`
    (This is not a state of the art algorithm)
    """
    primeFactors = []
    limit = int(math.sqrt(n)) + 1
    check = 2
    if n == 1: return [1]
    for check in range(2, limit):
        while n % check == 0:
            primeFactors.append(check)
            n /= check
    if n > 1:
        primeFactors.append(n)
    return primeFactors


def bruteForceV1(publicN, publicE, encrypted_msg):
    start = datetime.now()
    print ("bruteForceV1:\nFirst we factor n: {} to find it's two prime numbers . . .".format(publicN))
    (p, q) = primeFactorizationV1(publicN)
    print ("That took: {} find private (p, q) of n as {}".format(datetime.now() - start, (p, q)))
    public, uncoveredPrivateKey = generate_keypair(p, q, publicE)
    message = decrypt(uncoveredPrivateKey, encrypted_msg)
    print "Revealed message: {}".format(message)


def _backspace(message):
    return chr(8)*len(str(defaultP))


if __name__ == '__main__':
    print("RSA Encrypter/ Decrypter")

    # List some primes
    primes = []
    for x in range(29999900, 30000000):
        isPrime = is_prime(x)
        if isPrime:
            primes.append(str(x))
    print ("{} are prime!".format(", ".join(primes))

    defaultP = 29999947
    p = int(input("Enter a prime number (17, 19, 23, etc): {}".format(defaultP) + _backspace(defaultP) ) or defaultP)
    defaultQ = 29999999
    q = int(input("Enter another prime number (Not one you entered above): {}".format(defaultQ) + _backspace(defaultQ)) or defaultQ)
    print ("Generating your public/private keypairs now . . .")
    public, private = generate_keypair(p, q)
    print ("Your public key is ", public, " and your private key is ", private)
    defaultMessage = 'hello world'
    message = input("Enter a message to encrypt with your private key: {}".format(defaultMessage) + _backspace(defaultMessage)) or defaultMessage
    print ("Your message is: {}".format(message))
    encrypted_msg = encrypt(public, message)
    print ("Your encrypted message is: ")
    print (', '.join(map(str, encrypted_msg)))
    print ("Decrypting message with public key ", public, " . . .")
    print ("Your message is:")
    print (decrypt(private, encrypted_msg))

    print ("\nNow try to uncover the message with only information from the public key:")
    publicN = public[1]
    publicE = public[0]
    print ("n: {}".format(publicN))
    print ("e: {}".format(publicE))
    print ("and the encrypted message: {}".format(encrypted_msg))
    bruteForceV1(publicN, publicE, encrypted_msg)

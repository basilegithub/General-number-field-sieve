Thanks for looking at this project.

This is a project I started a few years ago, in 2019-2020.
I re-worked it recently to merge some versions, make it more readable, add comments, and get all the sources I used.


##### Introduction ######

This project is a Python implementation of the General number field sieve to factor integers.

Here are some detailed characteristics :
- Kleinjung 2006 algorithm to find good polynomials
- Double large prime variation used to collect relations for both algebraic and rational side
- Batch smoothness test and naive smoothness test available
- Union-find algorithm to detect cycles. Leveraging the graph structures (one central hub) to fasten the cycle finding when I know there is one.
- Parallel sieving allowed
- Parallel Polynomial search
- gaussian elimination, block Lanczos, Wiedemann algorithms available for the linear algebra step (no parallelization available for now)
- Lifting and Couveignes methods for computing the algebraic square root, no Montgomery algorithm available for now. (no parallelization for now)

##### Sources #####

Overall algorithm:
- "Prime numbers, a computational perspective" by Richard Crandall and Carl Pomerance (really good): https://link.springer.com/book/10.1007/0-387-28979-8
- "The Development of the Number Field Sieve" by many (really really good): https://link.springer.com/book/10.1007/BFb0091534
- "A beginner's guide to the general number field sieve" by Michael Case: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/NFSmadeeasy.pdf

Kleinjung polynomials search algorithm:
- "On polynomial selection for the number field sieve" by THORSTEN KLEINJUNG: https://www.ams.org/journals/mcom/2006-75-256/S0025-5718-06-01870-9/S0025-5718-06-01870-9.pdf

Polynomial metric Ep_score:
- "A new ranking function for polynomial selection in the number field sieve" by Nicolas David and Paul Zimmerman: https://inria.hal.science/hal-02151093v4/document

Double large prime
- "Factoring with two large primes" by Arjen K. Lenstra and Mark S. Manasse: https://scispace.com/pdf/factoring-with-two-large-primes-1lk9719aco.pdf

Batch smoothness test:
- "HOW TO FIND SMOOTH PARTS OF INTEGERS" by DANIEL J. BERNSTEIN: https://cr.yp.to/factorization/smoothparts-20040510.pdf

Gaussian elimination:
- Took the gaussian elimination code from https://github.com/skollmann/PyFactorise/blob/master/factorise.py#L39

Block Lanczos:
- "A Block Lanczos Algorithm for Finding Dependencies over GF(2)" by Peter L. Montgomery: https://scispace.com/pdf/a-block-lanczos-algorithm-for-finding-dependencies-over-gf-2-ezdu2qt0pp.pdf
- "A modified block Lanczos algorithm with fewer vectors" by Emmanuel Thomé (really good): https://eprint.iacr.org/2016/329.pdf

Wiedemann algorithm:
- "SOLVING HOMOGENEOUS LINEAR EQUATIONSOVER GF(2) VIA BLOCK WIEDEMANN ALGORITHM" by Don Coppersmith: https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1192970-7/S0025-5718-1994-1192970-7.pdf

Square root algorithms
- "Computing a square root for the number field sieve" by Jean-Marc Couveignes: https://www.math.u-bordeaux.fr/~jcouveig/publi/Cou94-2.pdf

##### running the algortihm #####

When running the main.py file, the required argument is --n the number you want to factor.

A suitable command looks like "python main.py --n [your_number]"

##### config parameters #####

I have defined a few parameters directly you can edit in the config file.

- batch_smooth_test: If True, the batch smoothness test from Bernstein will be used. Otherwise, the naive (trial division from the prime factor base) approach will be used. Batch smoothness test is generally faster

- gaussian_pivot: If True, gaussian elimination will ALWAYS be used in the linear algebra step. If you want to use block Lanczos of Wiedemann algorithms, this has to be set to False.

- lanczos: If True, the block Lanczos algorithm will be used. If False, the Wiedemann algorithm will be used. No matter its value, if gaussian_pivot is True then gaussian elimination will be performed.

- square_root_couveignes: If True, the Couveignes algorithm is run for the square root step. It requires many primes to be inert. If False, the lifting algorithm is used

- large_primes_constant: Define the constant that will be multiplied by the last prime in the factor base to obtain the bound for the single large primes. For the double large primes, the bound is this constant multiplied by the last prime in the factor base squared. Both the algebraic and rational sides have the same bounds.

- block_size: If Block Lanczos or Wiedemann algorithm are used, this sets the block size. In the Wiedemann algorithm, this is only used to compute many matrix-vector products at once, and not to compute the matrix generator of the matrix sequence you obtain.

- poly_search_nb_poly_coarse_eval: number of polynomials to be generated before doing precise ranking.

- poly_search_nb_poly_precise_eval: number of polynomials to be kept for precise evaluation.

- poly_search_prime_bound: maximum size of primes used in the polynomial generation

- poly_search_nb_roots: number of roots to use in the Kleinjung polynomial generation algorithm. It is the l parameter in the original paper.

- poly_search_multiplier: leading coefficient of generated polynomials is always a multiple of this multiplier.

- NB_CPU_POLY_SELECTION: Sets the number of cpu used for running polynomial search. One cpu is always kept as a "leader" that collect polynomial candidates from the workers.

- NB_CPU_SIEVE: Sets the number of cpu used for sieving. One cpu is always kept as a "leader" that collects relations from the sievers, and the rest of the cpus sieve, test for smoothness, and send their candidate relations to the leader.

##### General discussion #####

I tried overall to use as few exernal libraries as possible. While this makes the code way
	harder to write for me, and understand for you, I viewed it as a good challenge to
	write good code, and get to know hard and foreign concepts to me. It took a lot of
	work, but I got a bug free code, and got much more knowledgable on the subroutines
	thanks to this.

Here are some points I consider working on at some point:
- I still have some work to do to implement the state of the art polynomial search algorithms.
	Right now this generates decent polynomials, but it seems that for high degree, eg 5 or
	6) there are some useful techniques (rational approximation of root, translations and
	rotations).

- I have to implement the Montgomery algorithm for the algebraic square root computation.
	For the small numbers I have used my algorithm on, it is not a bottleneck yet.

- I have to optimize the Union-find algorithm for cycle detection, the version I have right now
	does not seem to be perfect, as I still do many operations to connect two distinct
	connected components.

- No matter how hard I tried, I am stuck on understanding the block Wiedemann algorithm. For
	now, the best I can do is the scalar one, with some optimizations. Namely, I use binary
	encoding of the blocks of vectors to compute matrix-vector product very efficiently.
	Then, for each scalar sequence, I compute its generator using scalar Berlekamp-Massey
	algorithm, and use it to update the minimal polynomial of the matrix. I consider
	working some way to compute the matrix generator of the matrix sequence through some
	algorithm (matrix Berlekamp-Massey or some other).

- I have to make the linear algebra step more parallel.

Here are the next steps:

- Some optimizations may be possible here and there, but I have already a big amount of work to
	get fast code. The next step on this topic is switching to C. Python is intrisically
	limited on the topic of speed, as it is compiled on the fly.
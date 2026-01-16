(defpackage :cma-es
  (:use :common-lisp :cffi)
  (:export #:make-optimizer
           #:ask
           #:tell
           #:run))

(in-package :cma-es)

                                        ; Loading the C Library

(cffi:define-foreign-library c-cmaes
  (:darwin "libcmaes.dylib")
  (:unix "libcmaes.so")
  (t (:default "libcmaes")))

(cffi:use-foreign-library c-cmaes)
  
                                        ; Foreign C Code

(defcfun ("cmaes_init" %cmaes-init) :pointer
  "Initializes the CMA-ES optimizer.
   Returns a pointer to a double buffer for allocating fitness scores."
  (evo :pointer)
  (dim :int)
  (xinit0 :pointer)
  (xstd0 :pointer)
  (seed :long)
  (popsize :int)
  (out :string))

(defcfun ("make_evo" %make-evo) :pointer
  "Makes a CMA-ES optimizer.
  Returns an opaque pointer to the cmaes_t struct.")


(defcfun ("cmaes_SamplePopulation" %ask) :pointer
  "Generates new candidate solutions."
  (evo :pointer))

(defcfun ("cmaes_UpdateDistribution" %tell) :void
  "Updates the distribution according to the fitnesses of the current
   population."
  (evo :pointer)
  (scores :pointer))

(defcfun ("cmaes_exit" destroy) :void
  "Frees the memory used by cmaes_t."
  (evo :pointer))
  
                                        ; Lisp Bindings

(defstruct cma-es-optimizer
  (evo-ptr) ; pointer to the cmaes_t struct
  (scores-ptr) ; pointer to a pre-allocated buffer for fitness scores
  (dim) ; the dimensionality of the problem
  (popsize)) ; the population size

(defun make-optimizer (dim pop-size &key
                                      (xinit0 (make-array dim :initial-element 0.5d0))
                                      (xstd0 (make-array dim :initial-element 0.5d0))
                                      (seed 0))
  "Makes a CMA-ES optimizer for tuning numerical constants."
  (cffi:with-foreign-objects ((means :double dim)
                              (stds :double dim))
    (loop for i from 0 below dim
          do (setf (cffi:mem-aref means :double i) (aref xinit0 i))
          do (setf (cffi:mem-aref stds :double i) (aref xstd0 i)))
    (let ((evo-ptr (%make-evo)))
      (make-cma-es-optimizer :evo-ptr evo-ptr
                             :scores-ptr (%cmaes-init evo-ptr dim means stds seed pop-size "none")
                             :dim dim
                             :popsize pop-size))))

(defun ask (optimizer)
  "Generate a new population of candidate solutions."
  (let* ((cmaes-t-struct (cma-es-optimizer-evo-ptr optimizer))
         (pop (%ask cmaes-t-struct))
         (pop-size (cma-es-optimizer-popsize optimizer))
         (dim (cma-es-optimizer-dim optimizer)))
    (loop for i from 0 below pop-size
          collect (loop for j from 0 below dim
                        collect (cffi:mem-aref (cffi:mem-aref pop :pointer i) :double j)))))

(defun tell (optimizer fitnesses)
  "Update the CMA-ES distribution according to the fitness of the population."
  (let ((pop-size (cma-es-optimizer-popsize optimizer))
        (scores (cma-es-optimizer-scores-ptr optimizer)))
    (loop for i from 0 below pop-size
          do (setf (cffi:mem-aref scores :double i) (coerce (elt fitnesses i) 'double-float)))
    (%tell (cma-es-optimizer-evo-ptr optimizer) scores)))

(defun run (fitness-fn
            dim 
            &key 
              (pop-size (+ 4 (floor (* 3 (log dim)))))
              (generations 100)
              (xinit0 (make-array dim :initial-element 0.5d0))
              (xstd0 (make-array dim :initial-element 0.5d0))
              (seed 0))
  (let ((opt (make-optimizer dim pop-size
                  :xinit0 xinit0
                  :xstd0 xstd0
                  :seed seed)))
    (loop repeat generations
          do (let* ((pop (ask opt))
                    (scores (mapcar fitness-fn pop)))
               (tell opt scores)))
    (let ((final-pop (ask opt)))
      (destroy (cma-es-optimizer-evo-ptr opt))
      (first (sort final-pop #'< :key fitness-fn)))))
                    
  

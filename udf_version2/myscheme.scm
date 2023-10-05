(define (make-new-rpvar name default type)
 (if (not (rp-var-object name))
 (rp-var-define name default type #f)))


(make-new-rpvar 'myudf/mass 0 'real)
(make-new-rpvar 'myudf/heat 0 'real)

(define gui-dialog-box
 (let ((dialog-box #f)
 (table)
 (myudf/box1)
 (myudf/box2)
 (myudf/mass)
 (myudf/heat)
 )

 (define (update-cb . args)
 (cx-set-real-entry myudf/mass (rpgetvar 'myudf/mass))
 (cx-set-real-entry myudf/heat (rpgetvar 'myudf/heat))
 )

 (define (apply-cb . args)
 (rpsetvar 'myudf/mass (cx-show-real-entry myudf/mass))
 (rpsetvar 'myudf/heat (cx-show-real-entry myudf/heat))
 (%run-udf-apply 2)
 ) 

 (lambda args
 (if (not dialog-box)
 (let ()
 (set! dialog-box (cx-create-panel "Mass and Heat Box" apply-cb update-cb))
 (set! table (cx-create-table dialog-box "" 'border #f 'below 0 'right-of 0))
 (set! myudf/box1 (cx-create-table table "mass box" 'row 0))
 (set! myudf/box2 (cx-create-table table "heat box" 'row 1))
 (set! myudf/mass (cx-create-real-entry myudf/box1 "For mass" 'row 0))
 (set! myudf/heat (cx-create-real-entry myudf/box2 "For heat" 'row 0))
 ) ;End Of Let Statement
 ) ;End Of If Statement
 ;Call To Open Dialog Box
 (cx-show-panel dialog-box)
 ) ;End Of Args Function
 ) ;End Of Let Statement
) ;End Of GUI-Dialog-Box Definition
(cx-add-menu "Control Mass and Heat" #f)
(cx-add-hitem "Control Mass and Heat" "Mass and Heat" #f #f #t #f)
(cx-add-item "Mass and Heat" "setup Box" #\U #f #t gui-dialog-box)
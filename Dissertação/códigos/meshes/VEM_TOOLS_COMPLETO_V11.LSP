;;; ================================================================
;;; VEM_TOOLS_COMPLETO_V11.LSP (Sem cálculo de n_int)
;;; ---------------------------------------------------------------
;;; Contem os comandos:
;;; - VEMEXPORT    : Exporta a malha (apenas nos fisicos: vertices e arestas).
;;; - VEMNODES     : Numera os nos no CAD e desenha um circulo no no.
;;; - VEMELEMENTS  : Numera os elementos (poligonos) no CAD.
;;; - VEMNODEERASE : Apaga textos e circulos gerados.
;;; ================================================================

;;; ================================================================
;;; COMANDO: VEMEXPORT
;;; ================================================================
(defun c:VEMEXPORT (/ 
                    dist2 find-node add-node get-vertices format-sci calc-edge-nodes calc-area-signed ensure-counterclockwise
                    tol k_vem dof_node ss i ent pts num_pts j pA pB edge_pts ep result nodes elements el_conn f filename conn_str idA idEP)

  (defun dist2 (p q)
    (+ (* (- (car p) (car q)) (- (car p) (car q)))
       (* (- (cadr p) (cadr q)) (- (cadr p) (cadr q))))
  )

  (defun find-node (p nodes tol / idx)
    (setq idx 0)
    (while (and (< idx (length nodes)) (> (dist2 p (nth idx nodes)) (* tol tol)))
      (setq idx (1+ idx))
    )
    (if (< idx (length nodes)) (1+ idx) nil)
  )

  (defun add-node (p n_type nodes tol / id)
    (setq id (find-node p nodes tol))
    (if id
      (list id nodes)
      (list (1+ (length nodes)) (append nodes (list (list (float (car p)) (float (cadr p)) n_type))))
    )
  )

  (defun calc-area-signed (pts / area n idx p1 p2)
    (setq area 0.0 n (length pts) idx 0)
    (while (< idx n)
      (setq p1 (nth idx pts))
      (if (= idx (1- n)) (setq p2 (nth 0 pts)) (setq p2 (nth (1+ idx) pts)))
      (setq area (+ area (- (* (car p1) (cadr p2)) (* (car p2) (cadr p1)))))
      (setq idx (1+ idx))
    )
    (/ area 2.0)
  )

  (defun ensure-counterclockwise (pts)
    (if (< (calc-area-signed pts) 0.0)
      (reverse pts)
      pts
    )
  )

  (defun get-vertices (ent tol / dados item xy raw_pts clean_pts prev p_first p_last)
    (setq dados (entget ent) raw_pts '())
    (foreach item dados
      (if (= (car item) 10)
        (progn
          (setq xy (cdr item))
          (setq raw_pts (append raw_pts (list (list (car xy) (cadr xy)))))
        )
      )
    )
    
    (setq clean_pts '())
    (if raw_pts
      (progn
        (setq prev (car raw_pts))
        (setq clean_pts (list prev))
        (foreach p (cdr raw_pts)
          (if (> (dist2 p prev) (* tol tol))
            (progn
              (setq clean_pts (append clean_pts (list p)))
              (setq prev p)
            )
          )
        )
      )
    )
    
    (if (> (length clean_pts) 2)
      (progn
        (setq p_first (car clean_pts) p_last (last clean_pts))
        (if (< (dist2 p_first p_last) (* tol tol))
          (setq clean_pts (reverse (cdr (reverse clean_pts))))
        )
      )
    )
    (ensure-counterclockwise clean_pts)
  )

  (defun format-sci (num / out)
    (setq out (rtos num 1 14))
    (vl-string-translate "e" "E" out)
  )

  (defun calc-edge-nodes (p1 p2 k / n pts_int dx dy p_int)
    (setq pts_int '())
    (if (> k 1)
      (progn
        (setq dx (- (car p2) (car p1)) dy (- (cadr p2) (cadr p1)) n 1)
        (while (< n k)
          (setq p_int (list (+ (car p1) (* dx (/ (float n) k))) (+ (cadr p1) (* dy (/ (float n) k)))))
          (setq pts_int (append pts_int (list p_int)))
          (setq n (1+ n))
        )
      )
    )
    pts_int
  )

  (setq k_vem (getint "\nOrdem VEM (k) <1>: "))
  (if (null k_vem) (setq k_vem 1))

  (setq dof_node (getint "\nGraus de liberdade por no <2>: "))
  (if (null dof_node) (setq dof_node 2))

  (setq tol (getreal "\nTolerancia para nos coincidentes <1e-8>: "))
  (if (null tol) (setq tol 1e-8))

  (princ "\nSelecione as LWPOLYLINEs dos elementos.")
  (setq ss (ssget '((0 . "LWPOLYLINE"))))
  (if (null ss) (setq ss (ssget "_X" '((0 . "LWPOLYLINE")))))

  (if (null ss)
    (princ "\nNenhuma LWPOLYLINE encontrada.")
    (progn
      (setq nodes '() elements '())
      
      ;; 1º Loop: Adiciona VÉRTICES (Tipo 1)
      (setq i 0)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3)
          (foreach p pts
            (setq result (add-node p 1 nodes tol))
            (setq nodes (cadr result))
          )
        )
        (setq i (1+ i))
      )

      ;; 2º Loop: Adiciona NÓS DE ARESTA (Tipo 2) se k > 1
      (if (> k_vem 1)
        (progn
          (setq i 0)
          (while (< i (sslength ss))
            (setq ent (ssname ss i))
            (setq pts (get-vertices ent tol))
            (if (>= (length pts) 3)
              (progn
                (setq num_pts (length pts) j 0)
                (while (< j num_pts)
                  (setq pA (nth j pts))
                  (if (= j (1- num_pts)) (setq pB (nth 0 pts)) (setq pB (nth (1+ j) pts)))
                  
                  (setq edge_pts (calc-edge-nodes pA pB k_vem))
                  (foreach ep edge_pts
                    (setq result (add-node ep 2 nodes tol))
                    (setq nodes (cadr result))
                  )
                  (setq j (1+ j))
                )
              )
            )
            (setq i (1+ i))
          )
        )
      )

      ;; 3º Loop: Conectividade dos Elementos
      (setq i 0)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3)
          (progn
            (setq el_conn '() num_pts (length pts) j 0)
            
            (while (< j num_pts)
              (setq pA (nth j pts))
              (if (= j (1- num_pts)) (setq pB (nth 0 pts)) (setq pB (nth (1+ j) pts)))
              
              (setq idA (find-node pA nodes tol))
              (setq el_conn (append el_conn (list idA)))
              
              (if (> k_vem 1)
                (progn
                  (setq edge_pts (calc-edge-nodes pA pB k_vem))
                  (foreach ep edge_pts
                    (setq idEP (find-node ep nodes tol))
                    (setq el_conn (append el_conn (list idEP)))
                  )
                )
              )
              (setq j (1+ j))
            )

            (setq elements (append elements (list el_conn)))
          )
        )
        (setq i (1+ i))
      )

      ;; Escreve o arquivo .dat
      (setq filename (getfiled "Salvar Malha VEM" "malha_vem" "dat" 1))
      (if filename
        (progn
          (setq f (open filename "w"))
          (write-line "! =====================================================================" f)
          (write-line "! 1. CABECALHO DA MALHA" f)
          (write-line "! [Nº de Nos Fisicos]  [Nº de Elementos]  [Ordem VEM (k)]  [DOFs por No]" f)
          (write-line "! =====================================================================" f)
          (write-line (strcat (itoa (length nodes)) "  " (itoa (length elements)) "  " (itoa k_vem) "  " (itoa dof_node)) f)
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 2. COORDENADAS E TIPOS DOS NOS (1=Vertice, 2=Aresta)" f)
          (write-line "! [X]  [Y]  [type_id]" f)
          (write-line "! =====================================================================" f)
          (setq i 1)
          (foreach n nodes
            (write-line (strcat "! No " (itoa i)) f)
            (write-line (strcat (format-sci (car n)) "  " (format-sci (cadr n)) "  " (itoa (caddr n))) f)
            (setq i (1+ i))
          )
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 3. CONECTIVIDADE DOS ELEMENTOS (Sentido Anti-Horario)" f)
          (write-line "! =====================================================================" f)
          (setq i 1)
          (foreach el elements
            (write-line (strcat "! Elemento " (itoa i) " (" (itoa (length el)) " nos de borda)") f)
            (write-line (itoa (length el)) f)
            (setq conn_str "")
            (foreach node_id el (setq conn_str (strcat conn_str (itoa node_id) "  ")))
            (write-line conn_str f)
            (setq i (1+ i))
          )
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 4. CONDICOES DE CONTORNO DE DIRICHLET (Deslocamentos Prescritos)" f)
          (write-line "! =====================================================================" f)
          (write-line "0" f)
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 5. MOLAS ELASTICAS CONCENTRADAS" f)
          (write-line "! =====================================================================" f)
          (write-line "0" f)

          (close f)
          (princ (strcat "\nExportacao concluida com sucesso! Arquivo salvo em: " filename))
        )
        (princ "\nExportacao cancelada.")
      )
    )
  )
  (princ)
)

;; Comandos auxiliares mantidos para suporte no AutoCAD
(defun c:VEMNODES (/ 
                   dist2 find-node add-node get-vertices calc-edge-nodes make-label make-circle calc-area-signed ensure-counterclockwise
                   k_vem tol h r_circle ss i ent pts num_pts j pA pB edge_pts ep result nodes id pt txt_pt)

  (defun dist2 (p q) (+ (* (- (car p) (car q)) (- (car p) (car q))) (* (- (cadr p) (cadr q)) (- (cadr p) (cadr q)))))
  (defun find-node (p nodes tol / idx)
    (setq idx 0)
    (while (and (< idx (length nodes)) (> (dist2 p (nth idx nodes)) (* tol tol))) (setq idx (1+ idx)))
    (if (< idx (length nodes)) (1+ idx) nil)
  )
  (defun add-node (p n_type nodes tol / id)
    (setq id (find-node p nodes tol))
    (if id (list id nodes) (list (1+ (length nodes)) (append nodes (list (list (float (car p)) (float (cadr p)) n_type)))))
  )
  (defun calc-area-signed (pts / area n idx p1 p2)
    (setq area 0.0 n (length pts) idx 0)
    (while (< idx n)
      (setq p1 (nth idx pts))
      (if (= idx (1- n)) (setq p2 (nth 0 pts)) (setq p2 (nth (1+ idx) pts)))
      (setq area (+ area (- (* (car p1) (cadr p2)) (* (car p2) (cadr p1)))))
      (setq idx (1+ idx))
    )
    (/ area 2.0)
  )
  (defun ensure-counterclockwise (pts)
    (if (< (calc-area-signed pts) 0.0) (reverse pts) pts)
  )
  (defun get-vertices (ent tol / dados item xy raw_pts clean_pts prev p_first p_last)
    (setq dados (entget ent) raw_pts '())
    (foreach item dados
      (if (= (car item) 10)
        (progn
          (setq xy (cdr item))
          (setq raw_pts (append raw_pts (list (list (car xy) (cadr xy)))))
        )
      )
    )
    (setq clean_pts '())
    (if raw_pts
      (progn
        (setq prev (car raw_pts))
        (setq clean_pts (list prev))
        (foreach p (cdr raw_pts)
          (if (> (dist2 p prev) (* tol tol))
            (progn
              (setq clean_pts (append clean_pts (list p)))
              (setq prev p)
            )
          )
        )
      )
    )
    (if (> (length clean_pts) 2)
      (progn
        (setq p_first (car clean_pts) p_last (last clean_pts))
        (if (< (dist2 p_first p_last) (* tol tol))
          (setq clean_pts (reverse (cdr (reverse clean_pts))))
        )
      )
    )
    (ensure-counterclockwise clean_pts)
  )
  (defun calc-edge-nodes (p1 p2 k / n pts_int dx dy p_int)
    (setq pts_int '())
    (if (> k 1)
      (progn (setq dx (- (car p2) (car p1)) dy (- (cadr p2) (cadr p1)) n 1) (while (< n k) (setq p_int (list (+ (car p1) (* dx (/ (float n) k))) (+ (cadr p1) (* dy (/ (float n) k))))) (setq pts_int (append pts_int (list p_int))) (setq n (1+ n)))))
    pts_int
  )
  (defun make-label (p txt h / ent)
    (entmakex (list '(0 . "TEXT") '(100 . "AcDbEntity") '(100 . "AcDbText") (cons 10 (list (car p) (cadr p) 0.0)) (cons 40 h) (cons 1 txt) (cons 7 (getvar "TEXTSTYLE")) '(72 . 1) '(73 . 2) (cons 11 (list (car p) (cadr p) 0.0))))
  )
  (defun make-circle (p r / ent)
    (entmakex (list '(0 . "CIRCLE") '(100 . "AcDbEntity") '(100 . "AcDbCircle") (cons 10 (list (car p) (cadr p) 0.0)) (cons 40 r) '(62 . 1)))
  )

  (setq k_vem (getint "\nOrdem VEM (k) <1>: "))
  (if (null k_vem) (setq k_vem 1))
  (setq tol (getreal "\nTolerancia para nos coincidentes <1e-8>: "))
  (if (null tol) (setq tol 1e-8))
  (setq h (getreal "\nAltura do texto <2.0>: "))
  (if (null h) (setq h 2.0))
  (setq r_circle (* h 0.25))

  (princ "\nSelecione as LWPOLYLINEs dos elementos para numerar NOS.")
  (setq ss (ssget '((0 . "LWPOLYLINE"))))
  (if (null ss) (setq ss (ssget "_X" '((0 . "LWPOLYLINE")))))

  (if (null ss)
    (princ "\nNenhuma LWPOLYLINE encontrada.")
    (progn
      (setq nodes '())
      (setq i 0)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3)
          (foreach p pts (setq result (add-node p 1 nodes tol)) (setq nodes (cadr result)))
        )
        (setq i (1+ i))
      )
      (if (> k_vem 1)
        (progn
          (setq i 0)
          (while (< i (sslength ss))
            (setq ent (ssname ss i))
            (setq pts (get-vertices ent tol))
            (if (>= (length pts) 3)
              (progn
                (setq num_pts (length pts) j 0)
                (while (< j num_pts)
                  (setq pA (nth j pts))
                  (if (= j (1- num_pts)) (setq pB (nth 0 pts)) (setq pB (nth (1+ j) pts)))
                  (setq edge_pts (calc-edge-nodes pA pB k_vem))
                  (foreach ep edge_pts (setq result (add-node ep 2 nodes tol)) (setq nodes (cadr result)))
                  (setq j (1+ j))
                )
              )
            )
            (setq i (1+ i))
          )
        )
      )
      (setq id 1)
      (foreach pt nodes
        (make-circle pt r_circle)
        (setq txt_pt (list (+ (car pt) (* h 0.4)) (+ (cadr pt) (* h 0.4))))
        (make-label txt_pt (itoa id) h)
        (setq id (1+ id))
      )
      (princ (strcat "\n" (itoa (length nodes)) " nos fisicos e circulos desenhados no CAD."))
    )
  )
  (princ)
)

(defun c:VEMELEMENTS (/ get-vertices dist2 make-label calc-centroid calc-area-signed ensure-counterclockwise tol h ss i ent pts centroid id)
  (defun dist2 (p q) (+ (* (- (car p) (car q)) (- (car p) (car q))) (* (- (cadr p) (cadr q)) (- (cadr p) (cadr q)))))
  (defun calc-area-signed (pts / area n idx p1 p2)
    (setq area 0.0 n (length pts) idx 0)
    (while (< idx n)
      (setq p1 (nth idx pts))
      (if (= idx (1- n)) (setq p2 (nth 0 pts)) (setq p2 (nth (1+ idx) pts)))
      (setq area (+ area (- (* (car p1) (cadr p2)) (* (car p2) (cadr p1)))))
      (setq idx (1+ idx))
    )
    (/ area 2.0)
  )
  (defun ensure-counterclockwise (pts) (if (< (calc-area-signed pts) 0.0) (reverse pts) pts))
  (defun get-vertices (ent tol / dados item xy raw_pts clean_pts prev p_first p_last)
    (setq dados (entget ent) raw_pts '())
    (foreach item dados
      (if (= (car item) 10)
        (progn
          (setq xy (cdr item))
          (setq raw_pts (append raw_pts (list (list (car xy) (cadr xy)))))
        )
      )
    )
    (setq clean_pts '())
    (if raw_pts
      (progn
        (setq prev (car raw_pts))
        (setq clean_pts (list prev))
        (foreach p (cdr raw_pts)
          (if (> (dist2 p prev) (* tol tol))
            (progn
              (setq clean_pts (append clean_pts (list p)))
              (setq prev p)
            )
          )
        )
      )
    )
    (if (> (length clean_pts) 2)
      (progn
        (setq p_first (car clean_pts) p_last (last clean_pts))
        (if (< (dist2 p_first p_last) (* tol tol))
          (setq clean_pts (reverse (cdr (reverse clean_pts))))
        )
      )
    )
    (ensure-counterclockwise clean_pts)
  )
  (defun calc-centroid (pts / sum_x sum_y n)
    (setq sum_x 0.0 sum_y 0.0 n (length pts))
    (foreach p pts (setq sum_x (+ sum_x (car p))) (setq sum_y (+ sum_y (cadr p))))
    (list (/ sum_x (float n)) (/ sum_y (float n)))
  )
  (defun make-label (p txt h / ent)
    (entmakex (list '(0 . "TEXT") '(100 . "AcDbEntity") '(100 . "AcDbText") (cons 10 (list (car p) (cadr p) 0.0)) (cons 40 h) (cons 1 txt) (cons 7 (getvar "TEXTSTYLE")) '(72 . 1) '(73 . 2) (cons 11 (list (car p) (cadr p) 0.0)) '(62 . 4)))
  )
  (setq tol (getreal "\nTolerancia para fechamento <1e-8>: "))
  (if (null tol) (setq tol 1e-8))
  (setq h (getreal "\nAltura do texto para os elementos <2.0>: "))
  (if (null h) (setq h 2.0))
  (princ "\nSelecione as LWPOLYLINEs dos elementos para numerar.")
  (setq ss (ssget '((0 . "LWPOLYLINE"))))
  (if (null ss) (setq ss (ssget "_X" '((0 . "LWPOLYLINE")))))
  (if (null ss)
    (princ "\nNenhuma LWPOLYLINE encontrada.")
    (progn
      (setq i 0 id 1)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3)
          (progn
            (setq centroid (calc-centroid pts))
            (make-label centroid (strcat "E" (itoa id)) h)
            (setq id (1+ id))
          )
        )
        (setq i (1+ i))
      )
      (princ (strcat "\n" (itoa (1- id)) " elementos numerados desenhados no CAD."))
    )
  )
  (princ)
)

(defun c:VEMNODEERASE (/ ss i)
  (princ "\nSelecione as numeracoes (textos) e circulos a apagar.")
  (setq ss (ssget '((0 . "TEXT,CIRCLE"))))
  (if ss
    (progn
      (setq i 0)
      (while (< i (sslength ss)) (entdel (ssname ss i)) (setq i (1+ i)))
      (princ (strcat "\n" (itoa (sslength ss)) " entidades apagadas."))
    )
    (princ "\nNenhuma entidade selecionada.")
  )
  (princ)
)

(princ "\nFerramentas VEM V11 carregadas.")
(princ)
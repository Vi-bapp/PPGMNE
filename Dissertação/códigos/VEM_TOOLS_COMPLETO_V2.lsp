;;; ================================================================
;;; VEM_TOOLS_COMPLETO_V2.LSP
;;; ---------------------------------------------------------------
;;; Contem os comandos VEMEXPORT e VEMNODES (e VEMNODEERASE).
;;; Ambos suportam ordem VEM (k) generica e DOFs configuraveis.
;;; CORRECAO: Nao exige mais a propriedade "Closed". Aceita qualquer
;;; polyline com >= 3 lados e limpa vertices sobrepostos (fechamento manual).
;;; ================================================================

;;; ================================================================
;;; COMANDO: VEMEXPORT
;;; ================================================================
(defun c:VEMEXPORT (/ 
                    dist2 find-node add-node get-vertices format-sci calc-edge-nodes
                    tol k_vem dof_node ss i ent pts num_pts j pA pB edge_pts ep result nodes elements el_conn f filename conn_str idA idEP)

  ;; --- Funcoes Internas ---
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

  (defun add-node (p nodes tol / id)
    (setq id (find-node p nodes tol))
    (if id
      (list id nodes)
      (list (1+ (length nodes)) (append nodes (list (list (float (car p)) (float (cadr p))))))
    )
  )

  ;; NOVA LOGICA: Pega os vertices, limpa o ultimo se for igual ao primeiro (fechamento manual)
  (defun get-vertices (ent tol / dados item xy pts_local p_first p_last)
    (setq dados (entget ent) pts_local '())
    (foreach item dados
      (if (= (car item) 10)
        (progn
          (setq xy (cdr item))
          (setq pts_local (append pts_local (list (list (car xy) (cadr xy)))))
        )
      )
    )
    (if (> (length pts_local) 2)
      (progn
        (setq p_first (car pts_local) p_last (last pts_local))
        (if (< (dist2 p_first p_last) (* tol tol))
          (setq pts_local (reverse (cdr (reverse pts_local)))) ; Remove o ultimo elemento se for repetido
        )
      )
    )
    pts_local
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

  ;; --- Logica Principal ---
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
      
      ;; PASSO 1: Coletar APENAS os vertices
      (setq i 0)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3) ;; Verifica se e um poligono valido (>= 3 lados)
          (foreach p pts
            (setq result (add-node p nodes tol))
            (setq nodes (cadr result))
          )
        )
        (setq i (1+ i))
      )

      ;; PASSO 2: Coletar APENAS os nos de aresta
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
                    (setq result (add-node ep nodes tol))
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

      ;; PASSO 3: Construir a Conectividade
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
              
              ;; Pega ID do vertice
              (setq idA (find-node pA nodes tol))
              (setq el_conn (append el_conn (list idA)))
              
              ;; Pega IDs dos nos da aresta
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

      ;; EXPORTACAO
      (setq filename (getfiled "Salvar Malha VEM" "malha_vem" "dat" 1))
      (if filename
        (progn
          (setq f (open filename "w"))
          (write-line "! =====================================================================" f)
          (write-line "! 1. CABECALHO DA MALHA" f)
          (write-line "! [Nº de Nos]  [Nº de Elementos]  [Ordem VEM (k)]  [DOFs por No]" f)
          (write-line "! =====================================================================" f)
          (write-line (strcat (itoa (length nodes)) "  " (itoa (length elements)) "  " (itoa k_vem) "  " (itoa dof_node)) f)
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 2. COORDENADAS E TIPOS DOS NOS (Primeiros = Vertices, Depois = Arestas)" f)
          (write-line "! [X]  [Y]  [type_id: 1=Vertice/Aresta]" f)
          (write-line "! =====================================================================" f)
          (setq i 1)
          (foreach n nodes
            (write-line (strcat "! No " (itoa i)) f)
            (write-line (strcat (format-sci (car n)) "  " (format-sci (cadr n)) "  1") f)
            (setq i (1+ i))
          )
          (write-line "" f)

          (write-line "! =====================================================================" f)
          (write-line "! 3. CONECTIVIDADE DOS ELEMENTOS (Poligonos/VEM)" f)
          (write-line "! =====================================================================" f)
          (setq i 1)
          (foreach el elements
            (write-line (strcat "! Elemento " (itoa i) " (" (itoa (length el)) " nos fisicos)") f)
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

;;; ================================================================
;;; COMANDO: VEMNODES (DESENHA NUMERACAO NO CAD)
;;; ================================================================
(defun c:VEMNODES (/ 
                   dist2 find-node add-node get-vertices calc-edge-nodes make-label
                   k_vem tol h ss i ent pts num_pts j pA pB edge_pts ep result nodes id pt)

  ;; --- Funcoes Internas ---
  (defun dist2 (p q) (+ (* (- (car p) (car q)) (- (car p) (car q))) (* (- (cadr p) (cadr q)) (- (cadr p) (cadr q)))))
  (defun find-node (p nodes tol / idx)
    (setq idx 0)
    (while (and (< idx (length nodes)) (> (dist2 p (nth idx nodes)) (* tol tol))) (setq idx (1+ idx)))
    (if (< idx (length nodes)) (1+ idx) nil)
  )
  (defun add-node (p nodes tol / id)
    (setq id (find-node p nodes tol))
    (if id (list id nodes) (list (1+ (length nodes)) (append nodes (list (list (float (car p)) (float (cadr p)))))))
  )
  
  (defun get-vertices (ent tol / dados item xy pts_local p_first p_last)
    (setq dados (entget ent) pts_local '())
    (foreach item dados
      (if (= (car item) 10)
        (progn
          (setq xy (cdr item))
          (setq pts_local (append pts_local (list (list (car xy) (cadr xy)))))
        )
      )
    )
    (if (> (length pts_local) 2)
      (progn
        (setq p_first (car pts_local) p_last (last pts_local))
        (if (< (dist2 p_first p_last) (* tol tol))
          (setq pts_local (reverse (cdr (reverse pts_local))))
        )
      )
    )
    pts_local
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

  ;; --- Logica Principal ---
  (setq k_vem (getint "\nOrdem VEM (k) <1>: "))
  (if (null k_vem) (setq k_vem 1))

  (setq tol (getreal "\nTolerancia para nos coincidentes <1e-8>: "))
  (if (null tol) (setq tol 1e-8))

  (setq h (getreal "\nAltura do texto <2.0>: "))
  (if (null h) (setq h 2.0))

  (princ "\nSelecione as LWPOLYLINEs dos elementos.")
  (setq ss (ssget '((0 . "LWPOLYLINE"))))
  (if (null ss) (setq ss (ssget "_X" '((0 . "LWPOLYLINE")))))

  (if (null ss)
    (princ "\nNenhuma LWPOLYLINE encontrada.")
    (progn
      (setq nodes '())
      
      ;; PASSO 1: Vertices
      (setq i 0)
      (while (< i (sslength ss))
        (setq ent (ssname ss i))
        (setq pts (get-vertices ent tol))
        (if (>= (length pts) 3)
          (foreach p pts (setq result (add-node p nodes tol)) (setq nodes (cadr result)))
        )
        (setq i (1+ i))
      )

      ;; PASSO 2: Arestas
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
                  (foreach ep edge_pts (setq result (add-node ep nodes tol)) (setq nodes (cadr result)))
                  (setq j (1+ j))
                )
              )
            )
            (setq i (1+ i))
          )
        )
      )

      ;; Desenha os textos
      (setq id 1)
      (foreach pt nodes
        (make-label pt (itoa id) h)
        (setq id (1+ id))
      )
      (princ (strcat "\n" (itoa (length nodes)) " nos numerados desenhados no CAD."))
    )
  )
  (princ)
)

;;; ================================================================
;;; COMANDO: VEMNODEERASE (APAGA OS TEXTOS GERADOS)
;;; ================================================================
(defun c:VEMNODEERASE (/ ss i)
  (princ "\nSelecione as numeracoes dos nos a apagar.")
  (setq ss (ssget '((0 . "TEXT"))))
  (if ss
    (progn
      (setq i 0)
      (while (< i (sslength ss))
        (entdel (ssname ss i))
        (setq i (1+ i))
      )
      (princ (strcat "\n" (itoa (sslength ss)) " textos apagados."))
    )
    (princ "\nNenhum texto selecionado.")
  )
  (princ)
)

(princ "\nFerramentas VEM carregadas. Comandos: VEMEXPORT, VEMNODES, VEMNODEERASE.")
(princ)

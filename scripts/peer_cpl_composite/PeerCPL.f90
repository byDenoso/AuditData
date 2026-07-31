    module PeerCPL
    use DarkEnergyInterface
    use Quintessence
    use results
    use classes
    implicit none

    private

    ! Composite dark-energy model with two independently evolved sectors:
    !   1. TEarlyQuintessence scalar field and its two perturbation variables;
    !   2. a late CPL component evolved with the CAMB PPF Gamma equation.
    !
    ! The scalar background integration includes the late CPL density in H(a),
    ! so this is not a post-processing sum of separately computed cosmologies.
    type, extends(TEarlyQuintessence) :: TPeerCPL
        real(dl) :: w = -1._dl
        real(dl) :: wa = 0._dl
        logical :: peer_enabled = .true.
        logical :: cpl_enabled = .true.
        real(dl) :: c_Gamma_ppf = 0.4_dl
        real(dl), private :: grhov0 = 0._dl
        integer, private :: peer_eqs = 0
        integer, private :: late_eqs = 0
        class(CAMBdata), pointer, private :: CompositeState => null()
    contains
        procedure :: Init => TPeerCPL_Init
        procedure :: BackgroundDensityAndPressure => TPeerCPL_BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy => TPeerCPL_PerturbedStressEnergy
        procedure :: PerturbationEvolve => TPeerCPL_PerturbationEvolve
        procedure :: diff_rhopi_Add_Term => TPeerCPL_diff_rhopi_Add_Term
        procedure :: Effective_w_wa => TPeerCPL_Effective_w_wa
        procedure :: EvolveBackground => TPeerCPL_EvolveBackground
        procedure, nopass :: PythonClass => TPeerCPL_PythonClass
        procedure, nopass :: SelfPointer => TPeerCPL_SelfPointer
        procedure, private :: LateBackground
        procedure, private :: PeerBackground
        procedure, private :: PeerPerturbations
    end type TPeerCPL

    public TPeerCPL

    contains

    function TPeerCPL_PythonClass()
    character(LEN=:), allocatable :: TPeerCPL_PythonClass

    TPeerCPL_PythonClass = 'PeerCPL'
    end function TPeerCPL_PythonClass


    subroutine TPeerCPL_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TPeerCPL), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TPeerCPL_SelfPointer


    subroutine TPeerCPL_Init(this, State)
    class(TPeerCPL), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    logical late_is_lambda

    select type(State)
    class is (CAMBdata)
        this%CompositeState => State
        this%grhov0 = State%grhov
    class default
        call MpiStop('PeerCPL requires CAMBdata state')
    end select

    this%peer_eqs = 0
    if (this%peer_enabled .and. this%fde_zc > 1e-12_dl) then
        ! The late component is supplied separately, so the scalar potential
        ! must not contain its own cosmological-constant floor.
        this%frac_lambda0 = 0._dl
        call this%TEarlyQuintessence%Init(State)
        this%peer_eqs = 2
    end if

    late_is_lambda = (.not. this%cpl_enabled) .or. &
        (abs(this%w + 1._dl) < 1e-12_dl .and. abs(this%wa) < 1e-12_dl)
    if (this%cpl_enabled .and. .not. late_is_lambda) then
        this%late_eqs = 1
    else
        this%late_eqs = 0
    end if

    this%c_Gamma_ppf = 0.4_dl
    this%num_perturb_equations = this%peer_eqs + this%late_eqs
    this%is_cosmological_constant = (this%peer_eqs == 0 .and. late_is_lambda)
    end subroutine TPeerCPL_Init


    subroutine LateBackground(this, a, grhov_t, wlate)
    class(TPeerCPL), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t, wlate
    real(dl) rel

    if (.not. this%cpl_enabled) then
        grhov_t = 0._dl
        wlate = -1._dl
        return
    end if

    wlate = this%w + this%wa * (1._dl - a)
    rel = a ** (1._dl - 3._dl * this%w - 3._dl * this%wa)
    if (this%wa /= 0._dl) rel = rel * exp(-3._dl * this%wa * (1._dl - a))
    grhov_t = this%grhov0 * rel / (a * a)
    end subroutine LateBackground


    subroutine PeerBackground(this, a, grhov_t, wpeer)
    class(TPeerCPL), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t, wpeer

    if (this%peer_eqs > 0) then
        call this%TQuintessence%BackgroundDensityAndPressure(this%grhov0, a, grhov_t, wpeer)
    else
        grhov_t = 0._dl
        wpeer = -1._dl
    end if
    end subroutine PeerBackground


    subroutine TPeerCPL_BackgroundDensityAndPressure(this, grhov, a, grhov_t, wtot)
    class(TPeerCPL), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: wtot
    real(dl) peer_rho, late_rho, peer_w, late_w, pressure

    call this%PeerBackground(a, peer_rho, peer_w)
    call this%LateBackground(a, late_rho, late_w)
    grhov_t = peer_rho + late_rho

    if (present(wtot)) then
        if (grhov_t > tiny(1._dl)) then
            pressure = peer_w * peer_rho + late_w * late_rho
            wtot = pressure / grhov_t
        else
            wtot = -1._dl
        end if
    end if
    end subroutine TPeerCPL_BackgroundDensityAndPressure


    subroutine TPeerCPL_EvolveBackground(this, num, a, y, yprime)
    ! Exact scalar-background evolution in the presence of the late CPL density.
    ! Variables are phi=y(1), a^2 phi'=y(2), matching TQuintessence.
    class(TPeerCPL), intent(in) :: this
    integer, intent(in) :: num
    real(dl), intent(in) :: a, y(num)
    real(dl), intent(out) :: yprime(num)
    real(dl) a2, phi, phidot, grho_peer_a4, grho_late_t, wlate, total, adot

    a2 = a * a
    phi = y(1)
    phidot = y(2) / a2
    grho_peer_a4 = a2 * (0.5_dl * phidot**2 + a2 * this%Vofphi(phi, 0))
    call this%LateBackground(a, grho_late_t, wlate)
    total = this%CompositeState%grho_no_de(a) + grho_peer_a4 + grho_late_t * a2
    adot = sqrt(total / 3._dl)
    yprime(1) = phidot / adot
    yprime(2) = -a2**2 * this%Vofphi(phi, 1) / adot
    end subroutine TPeerCPL_EvolveBackground


    subroutine PeerPerturbations(this, a, k, y, w_ix, dgrho_peer, dgq_peer)
    class(TPeerCPL), intent(in) :: this
    real(dl), intent(in) :: a, k
    real(dl), intent(in) :: y(:)
    integer, intent(in) :: w_ix
    real(dl), intent(out) :: dgrho_peer, dgq_peer
    real(dl) phi, phidot, clxq, vq

    if (this%peer_eqs == 0) then
        dgrho_peer = 0._dl
        dgq_peer = 0._dl
        return
    end if

    call this%ValsAta(a, phi, phidot)
    clxq = y(w_ix)
    vq = y(w_ix + 1)
    dgrho_peer = phidot * vq + clxq * a**2 * this%Vofphi(phi, 1)
    dgq_peer = k * phidot * clxq
    end subroutine PeerPerturbations


    subroutine TPeerCPL_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, wtot, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TPeerCPL), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, wtot, gpres_noDE, &
        etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) dgrho_peer, dgq_peer, dgrho_late, dgq_late
    real(dl) peer_rho, late_rho, peer_w, late_w, peer_pressure
    real(dl) Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
    real(dl) vT, grhoT, k2, dgq_other, dgrho_other, gpres_other
    integer late_ix

    call this%PeerPerturbations(a, k, ay, w_ix, dgrho_peer, dgq_peer)
    dgrho_late = 0._dl
    dgq_late = 0._dl

    if (this%late_eqs > 0) then
        call this%PeerBackground(a, peer_rho, peer_w)
        call this%LateBackground(a, late_rho, late_w)
        peer_pressure = peer_w * peer_rho
        late_ix = w_ix + this%peer_eqs
        k2 = k * k

        ! For the late PPF component, the physical scalar field belongs to
        ! the "other" sector in the PPF closure relations.
        dgrho_other = dgrho + dgrho_peer
        dgq_other = dgq + dgq_peer
        gpres_other = gpres_noDE + peer_pressure
        grhoT = grho - late_rho
        vT = dgq_other / (grhoT + gpres_other)
        Gamma = ay(late_ix)

        sigma = (etak + (dgrho_other + 3._dl * adotoa / k * dgq_other) / &
            (2._dl * k)) / kf1 - k * Gamma
        sigma = sigma / adotoa

        S_Gamma = late_rho * (1._dl + late_w) * (vT + sigma) * k / &
            (2._dl * adotoa * k2)
        ckH = this%c_Gamma_ppf * k / adotoa

        if (ckH * ckH > 1000._dl) then
            Gamma = 0._dl
            Gammadot = 0._dl
        else
            Gammadot = S_Gamma / (1._dl + ckH * ckH) - Gamma - ckH * ckH * Gamma
            Gammadot = Gammadot * adotoa
        end if
        ayprime(late_ix) = Gammadot

        Fa = 1._dl + 3._dl * (grhoT + gpres_other) / (2._dl * k2 * kf1)
        dgq_late = S_Gamma - Gammadot / adotoa - Gamma
        dgq_late = -dgq_late / Fa * 2._dl * k * adotoa + &
            vT * late_rho * (1._dl + late_w)
        dgrho_late = -2._dl * k2 * kf1 * Gamma - 3._dl / k * adotoa * dgq_late
    end if

    dgrhoe = dgrho_peer + dgrho_late
    dgqe = dgq_peer + dgq_late
    end subroutine TPeerCPL_PerturbedStressEnergy


    subroutine TPeerCPL_PerturbationEvolve(this, ayprime, wtot, w_ix, a, adotoa, k, z, y)
    class(TPeerCPL), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: wtot, a, adotoa, k, z, y(:)
    integer, intent(in) :: w_ix

    if (this%peer_eqs > 0) then
        call this%TQuintessence%PerturbationEvolve(ayprime, wtot, w_ix, a, adotoa, k, z, y)
    end if
    ! The PPF Gamma derivative is set in PerturbedStressEnergy, as in TDarkEnergyPPF.
    end subroutine TPeerCPL_PerturbationEvolve


    function TPeerCPL_diff_rhopi_Add_Term(this, a, dgrhoe, dgqe, grho, gpres, wtot, &
        grhok, adotoa, Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TPeerCPL), intent(in) :: this
    real(dl), intent(in) :: a, dgrhoe, dgqe, grho, gpres, wtot, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime(:), y(:)
    integer, intent(in) :: w_ix
    real(dl) ppiedot, hdotoh
    real(dl) dgrho_peer, dgq_peer, dgrho_late, dgq_late
    real(dl) late_rho, late_w
    integer late_ix

    if (this%late_eqs == 0) then
        ppiedot = 0._dl
        return
    end if

    call this%PeerPerturbations(a, k, y, w_ix, dgrho_peer, dgq_peer)
    dgrho_late = dgrhoe - dgrho_peer
    dgq_late = dgqe - dgq_peer
    call this%LateBackground(a, late_rho, late_w)
    late_ix = w_ix + this%peer_eqs

    hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / (6._dl * adotoa)
    ppiedot = 3._dl * dgrho_late + dgq_late * &
        (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
        late_rho * (1._dl + late_w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
        (yprime(late_ix) / adotoa - 2._dl * y(late_ix))
    ppiedot = ppiedot * adotoa / Kf1
    end function TPeerCPL_diff_rhopi_Add_Term


    subroutine TPeerCPL_Effective_w_wa(this, weff, waeff)
    class(TPeerCPL), intent(inout) :: this
    real(dl), intent(out) :: weff, waeff

    if (this%cpl_enabled) then
        weff = this%w
        waeff = this%wa
    else
        weff = -1._dl
        waeff = 0._dl
    end if
    end suboutine TPeerCPL_Effective_w_wa

    end module PeerCPLS

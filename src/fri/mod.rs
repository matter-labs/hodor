fn is_power_of_2(num: usize) -> bool {
    assert!(num > 0);

    let mut pow = 0;


    while (1 << (pow+1)) < num {
        pow += 1;
    }

    if (1 << (pow+1)) == num {
        return true;
    }

    false
}

mod reduction;
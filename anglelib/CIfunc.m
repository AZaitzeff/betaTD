function newCI=CIfunc(CI)
    newCI=double(CI>=.1)+(double(CI<.1).*double(CI>.03).*(CI-.03)/.07);
end
classdef Field < handle
    properties
        Wires
        eps
        N
        coordinates_center
        R
        Alpha
        tau
        fig
    end
    
    methods
        function obj = Field(Wires_data, eps)  % конструктор (аналог __init__)
            %FIELD Конструктор класса Field для моделирования электростатического поля.
            %   OBJ = FIELD(WIRES_DATA, EPS) инициализирует объект класса Field,
            %   формирует геометрию элементарных проводников с учетом зеркальных
            %   изображений относительно земли, рассчитывает матрицу потенциальных
            %   коэффициентов Alpha и вычисляет линейные плотности зарядов (tau).
            %
            %   Входные аргументы:
            %       WIRES_DATA - Массив структур с параметрами элементарных проводов.
            %                    Должен содержать обязательные поля:
            %                       .x  - координата центра провода по оси X (м)
            %                       .y  - координата центра провода по оси Y (м)
            %                       .U  - электрический потенциал провода (В)
            %                       .r0 - радиус элементарного проводника (м)
            %       EPS        - Относительная диэлектрическая проницаемость среды.
            %
            %   Выходные аргументы:
            %       OBJ        - Инициализированный handle-объект класса Field.
            %
            %   Вычисляемые свойства объекта:
            %       obj.N                  - количество реальных элементарных проводов
            %       obj.coordinates_center - матрица координат (2N x 2) для реальных 
            %                                и зеркальных проводников
            %       obj.Alpha              - матрица потенциальных коэффициентов (N x N)
            %       obj.tau                - вектор-столбец линейных зарядов (Кл/м)
            %
            %   Пример использования:
            %       field = Field(elementary_wires, 1);
          
            obj.Wires = Wires_data;
            obj.eps = eps;
        
            obj.N = length([obj.Wires.x]);
            obj.coordinates_center = [[obj.Wires.x]' [obj.Wires.y]';
                                  [obj.Wires.x]' -[obj.Wires.y]'];
            
            obj.R = [obj.Wires.r0];   % радиусы
            
            % Вычисляем матрицу расстояний между центрами
            D_centers = squareform(pdist(obj.coordinates_center));
            % Создание матрицы, где в строке i, столбце j стоит R_j
            R_matrix = repmat(obj.R, 2*obj.N, 1); % каждый столбец j содержит R_j
            % Формирование Alpha
            Distances = D_centers(:,1:obj.N) - R_matrix;
            % На диагонали должно быть R(i), а не 0 - R(i))
            Distances(1:2*obj.N+1:end) = obj.R;  % индексы диагональных элементов
            Alpha = [log(Distances((obj.N+1):end, :) ./ Distances(1:obj.N, :))]';
            
            obj.Alpha = 1/(2*pi*obj.eps*8.8541878128e-12)*Alpha;
            obj.tau = obj.Alpha\[obj.Wires.U]';
        end
        
        function E = E_calculation(obj, x, y)  % метод
            %E_CALCULATION Расчет модуля напряженности электростатического поля в точке.
            %   E = E_CALCULATION(OBJ, X, Y) вычисляет модуль вектора напряженности
            %   электрического поля E (В/м) в точке с координатами (X, Y) методом
            %   заркальных отображений.
            %
            %   Граничные и физические условия:
            %       - Если точка находится на земле или под ней (Y < 0), метод 
            %         возвращает E = 0.
            %       - Если точка попадает внутрь или на границу любого из элементарных
            %         проводов, метод возвращает E = 0.
            %
            %   Входные аргументы:
            %       OBJ - Объект класса Field.
            %       X   - Координата запрашиваемой точки по оси X (м).
            %       Y   - Координата запрашиваемой точки по оси Y (м).
            %
            %   Выходные аргументы:
            %       E   - Модуль напряженности электрического поля (В/м).
            %
            %   Пример использования:
            %       E_val = field.E_calculation(10, 5);

            % Проверка под землей напряженность равна 0
            if y < 0
                E = 0;
                return;
            end

            % внутри любого провода E = 0
            % Находим расстояния от точки (x, y) до всех центров 
            if any(hypot(obj.coordinates_center(1:obj.N, 1) - x, obj.coordinates_center(1:obj.N, 2) - y) <= obj.R(:))
                E = 0;
                return;
            end
            
            dx = obj.coordinates_center(1:obj.N, 1) - x;
            dy_real = obj.coordinates_center(1:obj.N, 2) - y;
            dy_image = obj.coordinates_center((obj.N+1):end, 2) - y;
            
            r2_real = dx.^2 + dy_real.^2;
            r2_image = dx.^2 + dy_image.^2;
            
            Ex = sum(obj.tau .* dx .* (1./r2_real - 1./r2_image)) * 1/(2*pi*obj.eps*8.8541878128e-12);
            Ey = sum(obj.tau .* (dy_real./r2_real - dy_image./r2_image)) * 1/(2*pi*obj.eps*8.8541878128e-12);

            E = sqrt(Ex^2 + Ey^2); % Искомое значение напряженности в точке А
        end

        function Phi = Phi_calculation(obj, x, y)
            %PHI_CALCULATION Расчет электрического потенциала в точке пространства.
            %   PHI = PHI_CALCULATION(OBJ, X, Y) вычисляет значение электрического
            %   потенциала phi (В) в точке с координатами (X, Y) методом
            %   зеркальных отображений.
            %
            %   Граничные и физические условия:
            %       - Если точка находится на земле или под ней (Y <= 0), метод 
            %         возвращает Phi = 0 (поверхность земли принята за ноль).
            %       - Если точка попадает внутрь или на границу любого из элементарных
            %         проводов, метод возвращает точный заданный потенциал этого провода U.
            %
            %   Входные аргументы:
            %       OBJ - Объект класса Field.
            %       X   - Координата запрашиваемой точки по оси X (м).
            %       Y   - Координата запрашиваемой точки по оси Y (м).
            %
            %   Выходные аргументы:
            %       PHI - Значение электрического потенциала (В).
            %
            %   Пример использования:
            %       phi_val = field.Phi_calculation(10, 5);


            %Проверка: на земле и под ней потенциал равен 0
            if y <= 0
                Phi = 0;
                return;
            end

            % внутри любого провода Phi = U_провода
            % Находим логический вектор: где точка попала внутрь провода
            inside_indx = hypot(obj.coordinates_center(1:obj.N, 1) - x, obj.coordinates_center(1:obj.N, 2) - y) <= obj.R(:);
            
            if any(inside_indx)
                Phi = obj.Wires(find(inside_indx, 1)).U;
                return;
            end


            dx = obj.coordinates_center(1:obj.N, 1) - x;
            dy_real = obj.coordinates_center(1:obj.N, 2) - y;
            dy_image = obj.coordinates_center((obj.N+1):end, 2) - y;
            
            r2_real = dx.^2 + dy_real.^2;
            r2_image = dx.^2 + dy_image.^2;
            
            Phi = sum(obj.tau .* log(sqrt(r2_image ./ r2_real))) * 1/(2*pi*obj.eps*8.8541878128e-12);
        end

        function plot_wires(obj)
            %PLOT_WIRES Визуализация геометрии проводов и заземляющей плоскости.
            %   PLOT_WIRES(OBJ) создает графическое окно, настраивает координатную
            %   сетку, отображает заземляющую плоскость (Y = 0) и все элементарные
            %   проводники в виде окружностей с указанием их потенциалов.
            %
            %
            %   Входные аргументы:
            %       OBJ - Объект класса Field (handle).
            %
            %   Пример использования:
            %       field.plot_wires(); 
            
            % Создание графического окна
            obj.fig = figure('Name', 'Координатная сетка', 'Position', [2000 100 1200 800]);
            
            % Расчет пределов осей координат
            x_min = min(obj.coordinates_center(:, 1)) - 5;
            x_max = max(obj.coordinates_center(:, 1)) + 5;
            y_min = -1;
            y_max = max(obj.coordinates_center(:, 2)) + 5;

            % Активация режима наложения и отображения сетки
            hold on;
            grid on;
            axis([x_min x_max y_min y_max]);
            xlabel('X');
            ylabel('Y');
            title('Кликните мышью в любом месте для вычисления f(x,y)');
            
            % Отрисовка линии поверхности земли
            ground_x = linspace(x_min, x_max);
            ground_y = zeros(1, length(ground_x));
            plot(ground_x, ground_y, 'vk', 'LineWidth', 0.5, 'MarkerFaceColor', 'k');
        
            % Настройка параметров отображения осей и фона фигуры
            set(gca, 'FontSize', 12);
            set(gcf, 'Color', 'white');
            
            % Настройка стиля и прозрачности линий сетки
            set(gca, 'GridLineStyle', '--');
            set(gca, 'GridAlpha', 0.6);

            % Отрисовка геометрии каждого элементарного проводника
            for i = 1:obj.N
                
                % Генерация точек окружности проводника
                theta = linspace(0, 2*pi, 100);
                x_circle = obj.coordinates_center(i, 1) + obj.R(i) * cos(theta);
                y_circle = obj.coordinates_center(i, 2) + obj.R(i) * sin(theta);
                
                % Выделение тела проводника цветом и нанесение маркера центра
                fill(x_circle, y_circle, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
                plot(obj.coordinates_center(i, 1), obj.coordinates_center(i, 2), 'b+', 'MarkerSize', 8, 'LineWidth', 1.5);

                % Вывод текстовой метки потенциала над выбранными проводниками
                if obj.Wires(i).show_label
                    text(obj.coordinates_center(i, 1), obj.coordinates_center(i, 2) + obj.R(i) + 0.1, sprintf('\\phi = %.2f', obj.Wires(i).U), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', ...
                        'FontSize', 9, ...
                        'Color', 'black', ...
                        'FontWeight', 'bold');
                end
            end
            
            % Нанесение финальных подписей осей и заголовка графика
            xlabel('X, м');
            ylabel('Y, м');
            title('Координатная сетка и система проводов');

            % Добавление аннотации с единицами измерения физических величин
            annotation(obj.fig, 'textbox', [0.74, 0.80, 0.15, 0.10], ...
                'String', sprintf('Единицы измерения:\nE: В/м\n\\phi: В'), ...
                'FontSize', 10, ...
                'FontWeight', 'normal', ...
                'BackgroundColor', [1, 1, 1], ... 
                'EdgeColor', [0.7, 0.7, 0.7], ...  
                'LineWidth', 0.5, ...
                'FitBoxToText', 'on', ...          
                'Margin', 6);                      

            % Деактивация режима наложения графических объектов
            hold off;
        end

        function plot_E_phi_at_point(obj, point_coordinates, plot_var)
            %PLOT_E_PHI_AT_POINT Визуализация расчетной точки и параметров поля на графике.
            %   PLOT_E_PHI_AT_POINT(OBJ, POINT_COORDINATES, PLOT_VAR) активирует
            %   текущую фигуру OBJ.FIG, наносит маркер расчетной точки с заданными
            %   координатами и добавляет текстовые подписи значений параметров поля.
            %
            %   Входные аргументы:
            %       OBJ               - Объект класса Field (handle).
            %       POINT_COORDINATES - Вектор-строка из 2 элементов [X, Y] (м).
            %       PLOT_VAR          - Направление расчета параметров поля:
            %                           "E"    - вывод только напряженности (В/м),
            %                           "phi"  - вывод только потенциала (В),
            %                           "both" - вывод обоих параметров.
            %
            %   Пример использования:
            %       field.plot_E_phi_at_point([1.5, 3.0], "both");

            arguments
                obj
                point_coordinates (1,2) {mustBeNumeric}
                plot_var {mustBeMember(plot_var, {'E', 'phi', 'both'})}
            end

            % Проверка существования графического окна
            if isempty(obj.fig) || ~isvalid(obj.fig)
                error('График не найден. Сначала вызовите метод plot_wires().');
            end

            % Активация сохраненной фигуры и режима наложения объектов
            figure(obj.fig);
            hold on;

            % Нанесение маркера расчетной точки на график
            plot(point_coordinates(1), point_coordinates(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            
            % Вывод текстовых меток параметров поля в зависимости от выбора пользователя
            if plot_var == "E"
                text(point_coordinates(1), point_coordinates(2), sprintf('  {\\itE}=%.3f', obj.E_calculation(point_coordinates(1), point_coordinates(2))), ...
                'FontSize', 9, 'Color', 'blue', 'FontWeight', 'bold');
            elseif plot_var == "phi"
                text(point_coordinates(1), point_coordinates(2), sprintf('  \\phi=%.3f', obj.Phi_calculation(point_coordinates(1), point_coordinates(2))), ...
                'FontSize', 9, 'Color', 'blue', 'FontWeight', 'bold');
            else
                text(point_coordinates(1), point_coordinates(2), sprintf('  {\\itE}=%.3f\n  \\phi=%.3f ', obj.E_calculation(point_coordinates(1), point_coordinates(2)), obj.Phi_calculation(point_coordinates(1), point_coordinates(2))), ...
                'FontSize', 9, 'Color', 'blue', 'FontWeight', 'bold');
            end

            % Деактивация режима наложения графических объектов
            hold off;
        end

        function point_coordinates = get_point_from_fig(obj)
            %GET_POINT_FROM_FIG Интерактивное получение координат точки с графика.
            %   POINT_COORDINATES = GET_POINT_FROM_FIG(OBJ) приостанавливает работу
            %   программы, активирует графическое окно OBJ.FIG и переводит курсор
            %   мыши в режим ожидания клика пользователя. Возвращает координаты 
            %   выбранной точки.
            %
            %   Входные аргументы:
            %       OBJ - Объект класса Field (handle).
            %
            %   Выходные аргументы:
            %       POINT_COORDINATES - Вектор-строка из 2 элементов [X, Y] (м).
            %                           Содержит метрические координаты клика
            %                           в масштабе текущих осей графика.

            
            % Проверка существования графического окна
            if isempty(obj.fig) || ~isvalid(obj.fig)
                error('График не найден. Сначала вызовите метод plot_wires().');
            end
            
            % Активация сохраненной фигуры перед выполнением клика мыши
            figure(obj.fig);
            
            % Считывание координат одного интерактивного клика пользователя
            [x, y] = ginput(1);
            
            % Формирование вектора-строки для совместимости с методом plot_E_phi_at_point
            point_coordinates = [x, y];
        end

    end
end